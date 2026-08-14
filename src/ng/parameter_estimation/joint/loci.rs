//! Which loci every sample keeps raw evidence at.
//!
//! The samples are walked separately — on different machines, at different times, with no
//! sample able to see what any other chose — so the set cannot be negotiated or handed
//! round. **This file is the rule that lets each sample arrive at the identical set on its
//! own**, for the two kinds of locus the two paths care about: ordinary positions for the
//! SNP/indel path, and repeat tracts for the STR path.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_loci.md`. Types:
//! `doc/devel/ng/arch/parameter_prepass_joint_loci.md`.
//!
//! Two rules, and only the first lives here:
//!
//! - **generic** — keep position `p` when `hash(contig, p, seed) < threshold`, the
//!   threshold set from the selectable length so that about `target` positions survive;
//! - **STR** — the `cap` lowest hashes within each *(period, repeat count)* stratum, which
//!   is [`RepeatCatalog::sample_loci_per_stratum`](crate::ng::repeat_catalog::RepeatCatalog::sample_loci_per_stratum)
//!   and ships with the catalog. This unit states the policy and calls it.
//!
//! **Nothing here opens an alignment.** The selection is a pure function of the reference,
//! the selectable regions, the catalog, a seed and two caps — which is what makes it
//! testable by running it twice and comparing.

use std::hash::{Hash, Hasher};

use md5::{Digest, Md5};
use xxhash_rust::xxh3::Xxh3;

use crate::ng::reference_info::{ContigInfo, ReferenceBasesObserver, ReferenceInfo};
use crate::ng::region_typing::RegionKind;
use crate::ng::repeat_catalog::{
    ReadScope, RepeatCatalog, RepeatCatalogError, StrRepeatCriteria, StratumCounts, StratumSample,
};
use crate::ng::tandem_repeat::ScanParams;
use crate::ng::types::{Bp, ContigId, GenomePosition, GenomeRegion, Position};

// ---------------------------------------------------------------------
// Where the selection may look
// ---------------------------------------------------------------------

/// The stretches of reference a run may keep positions from, and their total length.
///
/// **One value carries both, deliberately.** The threshold that yields about `target`
/// positions is `2^64 · target / length`, and the sweep that applies it walks these same
/// stretches. Two statements of "how much genome is in play" — a contig table for the
/// arithmetic and a region list for the sweep — is the failure
/// `parameter_prepass_joint_loci.md` §7.1 names as the one most likely to be written by
/// reflex: under a `--regions` BED the count then comes out of the genome's length and
/// every estimate is scaled by the ratio. Here it cannot, because [`total_length`] is a
/// property of the same object the sweep reads.
///
/// [`total_length`]: SelectableRegions::total_length
///
/// **Two things narrow it, and they are different in kind.** The analysed regions are the
/// run's own choice — with no `--regions` BED, every contig whole. Reference `N` runs are
/// not a choice: no read can align there, so a position inside one is kept, never
/// observed, and enters every per-locus rate as a denominator with no numerator. Both are
/// resolved before this type is built, so the rule below never has to know which is which.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectableRegions {
    /// Ascending by `(contig, start)`, non-empty, and non-overlapping — the invariant
    /// [`new`](SelectableRegions::new) establishes and everything below relies on for the
    /// kept positions to come out sorted and duplicate-free.
    spans: Vec<GenomeRegion>,
    total_length: Bp,
}

impl SelectableRegions {
    /// Sort, check and total the spans.
    ///
    /// # Errors
    ///
    /// [`SelectionError::OverlappingRegions`] when two spans on one contig overlap after
    /// sorting, and [`SelectionError::InvertedRegion`] when a span ends before it starts.
    /// Both would let one position be offered twice, and a duplicated kept position is a
    /// silent double-weight in every rate fitted from the set.
    pub fn new(mut spans: Vec<GenomeRegion>) -> Result<Self, SelectionError> {
        spans.sort_unstable_by_key(|r| (r.contig.get(), r.start.get(), r.end.get()));
        let mut total: u64 = 0;
        for (i, span) in spans.iter().enumerate() {
            if span.end.get() < span.start.get() {
                return Err(SelectionError::InvertedRegion { region: *span });
            }
            if i > 0 {
                let previous = spans[i - 1];
                if previous.contig == span.contig && span.start.get() <= previous.end.get() {
                    return Err(SelectionError::OverlappingRegions {
                        first: previous,
                        second: *span,
                    });
                }
            }
            total += span.end.get() - span.start.get() + 1;
        }
        Ok(Self {
            spans,
            total_length: Bp(total),
        })
    }

    /// Every stretch, ascending by `(contig, start)`.
    pub fn spans(&self) -> &[GenomeRegion] {
        &self.spans
    }

    /// How many bases the selection may draw from — **the denominator of
    /// [`threshold_for`]**, and never a contig table's total.
    pub fn total_length(&self) -> Bp {
        self.total_length
    }

    /// The stretches both sets hold.
    ///
    /// **The one operation that narrows a domain**, and it returns the same type so that
    /// the narrowed length travels with the narrowed spans. Intersecting two sorted,
    /// disjoint sets cannot produce an overlap or an inversion, so the result is built
    /// directly rather than through [`new`](Self::new) — there is nothing for that check to
    /// find, and the total is summed as the pieces are cut.
    pub fn intersect(&self, other: &Self) -> Self {
        let mut spans = Vec::new();
        let mut total: u64 = 0;
        let (mut mine, mut theirs) = (0, 0);
        while mine < self.spans.len() && theirs < other.spans.len() {
            let (left, right) = (self.spans[mine], other.spans[theirs]);
            if left.contig != right.contig {
                // Both lists are sorted by contig first, so advance whichever is behind.
                if left.contig.get() < right.contig.get() {
                    mine += 1;
                } else {
                    theirs += 1;
                }
                continue;
            }
            let start = left.start.get().max(right.start.get());
            let end = left.end.get().min(right.end.get());
            if start <= end {
                spans.push(GenomeRegion {
                    contig: left.contig,
                    start: Position(start),
                    end: Position(end),
                });
                total += end - start + 1;
            }
            // Retire whichever span ends first; the other may still meet the next one.
            if left.end.get() < right.end.get() {
                mine += 1;
            } else {
                theirs += 1;
            }
        }
        Self {
            spans,
            total_length: Bp(total),
        }
    }
}

/// The stretches of a reference that are actually `A`, `C`, `G` or `T`, collected as the
/// reference streams past.
///
/// **Why the selection needs this at all.** The generic rule is a hash of a coordinate, so
/// it happily keeps positions inside an assembly gap. Such a position is never covered by
/// a read in any sample, so it contributes no evidence and yet sits in the denominator of
/// every rate derived from the kept set — where the posterior at a locus with no reads is
/// the prior, and the "rate" it contributes is the model's own prediction rather than a
/// measurement.
///
/// **Measured, and it is a property of the reference rather than a rounding term**: with a
/// two-million-position budget, tomato SL4.00 puts **135** kept positions inside `N` and
/// GRCh38 with the hs38d1 decoys puts **106,423** — 5.3% of the whole budget, because 5.3%
/// of that assembly is `N` (`examples/ng_joint_loci_probe.rs`).
///
/// **It costs no read of its own.** This is a [`ReferenceBasesObserver`], so the mask comes
/// off the same forward pass that computes the reference's digests — the seam the repeat
/// catalog uses for the same reason.
#[derive(Debug, Default)]
pub struct UnambiguousRuns {
    /// One entry per contig, in reference order: maximal runs of unambiguous base, 1-based
    /// and inclusive.
    runs: Vec<Vec<(u64, u64)>>,
    ambiguous_bases: u64,
    /// Where the current contig's next base sits, 1-based.
    next_position: u64,
    /// The run being extended, if any.
    open_run: Option<(u64, u64)>,
}

impl UnambiguousRuns {
    /// How many bases were not `A`, `C`, `G` or `T`.
    pub fn ambiguous_bases(&self) -> u64 {
        self.ambiguous_bases
    }

    /// The runs as a selection domain, ready for [`threshold_for`] and
    /// [`select_generic_positions`].
    ///
    /// # Errors
    ///
    /// Only what [`SelectableRegions::new`] refuses, which maximal runs cannot trip — the
    /// call is fallible because the constructor is, not because this can fail.
    pub fn into_selectable(self) -> Result<SelectableRegions, SelectionError> {
        SelectableRegions::new(
            self.runs
                .iter()
                .enumerate()
                .flat_map(|(i, contig_runs)| {
                    contig_runs.iter().map(move |&(start, end)| GenomeRegion {
                        contig: ContigId(i as u32),
                        start: Position(start),
                        end: Position(end),
                    })
                })
                .collect(),
        )
    }

    fn close_run(&mut self) {
        if let Some(run) = self.open_run.take() {
            self.runs
                .last_mut()
                .expect("a contig is open whenever a run is")
                .push(run);
        }
    }
}

impl ReferenceBasesObserver for UnambiguousRuns {
    fn contig_started(&mut self, _name: &str, _index: usize) {
        self.runs.push(Vec::new());
        self.next_position = 1;
        self.open_run = None;
    }

    fn bases(&mut self, upper: &[u8]) {
        for base in upper {
            let position = self.next_position;
            self.next_position += 1;
            if matches!(base, b'A' | b'C' | b'G' | b'T') {
                match &mut self.open_run {
                    Some(run) if run.1 + 1 == position => run.1 = position,
                    _ => {
                        self.close_run();
                        self.open_run = Some((position, position));
                    }
                }
            } else {
                self.ambiguous_bases += 1;
                self.close_run();
            }
        }
    }

    fn contig_finished(&mut self, _info: &ContigInfo) {
        self.close_run();
    }
}

// ---------------------------------------------------------------------
// The generic rule
// ---------------------------------------------------------------------

/// The value a position is selected by: a hash of **where it is**, not of what is there.
///
/// Position and not content, so that the selection cannot come to depend on the data; and
/// a fixed hash rather than a random number generator, so the same reference, region set
/// and seed keep the identical positions on every machine and in every run.
///
/// **Keyed by the contig's name, and by the same construction the STR sampler uses**
/// (`repeat_catalog::strata::hash_locus`), so both halves of the selection rest on one
/// uniformity assumption rather than two. A name rather than an index because an index is
/// a property of the file's contig order: two references that hold the same contigs in a
/// different order would select different positions, and the reference digest would say
/// so — but only after the fact.
pub fn hash_position(contig_name: &str, position: u64, seed: u64) -> u64 {
    ContigHasher::new(contig_name, seed).hash_position(position)
}

/// [`hash_position`] with the contig's name already absorbed.
///
/// One contig's name is hashed once and the per-position work is a clone and one `u64`.
/// **The value is identical to [`hash_position`]'s** — the same bytes reach the same
/// hasher in the same order — which is what the `hoisting_the_contig_name_changes_no_hash`
/// test pins, because the whole selection rests on those two agreeing.
#[derive(Clone)]
pub struct ContigHasher {
    base: Xxh3,
}

impl ContigHasher {
    pub fn new(contig_name: &str, seed: u64) -> Self {
        let mut base = Xxh3::with_seed(seed);
        contig_name.hash(&mut base);
        Self { base }
    }

    #[inline]
    pub fn hash_position(&self, position: u64) -> u64 {
        let mut hasher = self.base.clone();
        position.hash(&mut hasher);
        hasher.finish()
    }
}

/// The threshold that yields about `target` positions out of `selectable_length`.
///
/// The hash takes `2^64` values and spreads uniformly across them, so keeping the fraction
/// `target / selectable_length` of them means keeping those below
/// `2^64 · target / selectable_length`. **Computed in 128 bits**: `2^64 · target`
/// overflows a `u64` for any target above zero. Saturates to "keep everything" when
/// `target >= selectable_length`, which is the right answer rather than a case to guard.
///
/// The realised count is binomial around `target`, so about `±√target` — roughly ±1,400 at
/// two million. Nothing downstream needs the count to be exact; what it needs is that
/// every sample realises the *same* count, which follows from the rule being a function of
/// position alone.
pub fn threshold_for(target: u64, selectable_length: Bp) -> u64 {
    let length = selectable_length.get();
    if length == 0 || target >= length {
        return u64::MAX;
    }
    let scaled = (u128::from(target) << 64) / u128::from(length);
    // `target < length` gives a quotient strictly below `2^64`, so the cast is exact; the
    // clamp is the belt to that argument's braces.
    u64::try_from(scaled).unwrap_or(u64::MAX)
}

/// Keep this position? — the whole generic rule, on one position.
#[inline]
pub fn keeps_position(hasher: &ContigHasher, position: Position, threshold: u64) -> bool {
    hasher.hash_position(position.get()) < threshold
}

/// Every position the generic rule keeps, ascending in genome order.
///
/// `contig_names` is indexed by [`ContigId`], as the reference's contig table is.
///
/// **Order-independent by construction**: the rule reads one position and the seed, so a
/// region-sharded sweep and a single-threaded one keep the same positions and merging is
/// concatenation in genome order.
///
/// # Errors
///
/// [`SelectionError::UnknownContig`] when a span names a contig the reference table does
/// not hold — a silent skip there is a whole chromosome missing from the selection, and
/// every sample missing it identically, so nothing downstream could notice.
pub fn select_generic_positions(
    regions: &SelectableRegions,
    contig_names: &[String],
    seed: u64,
    threshold: u64,
) -> Result<Vec<GenomePosition>, SelectionError> {
    let mut kept = Vec::new();
    let mut current: Option<(ContigId, ContigHasher)> = None;
    for span in regions.spans() {
        let hasher = match &current {
            Some((contig, hasher)) if *contig == span.contig => hasher,
            _ => {
                let name = contig_names.get(span.contig.get() as usize).ok_or(
                    SelectionError::UnknownContig {
                        contig: span.contig,
                    },
                )?;
                current = Some((span.contig, ContigHasher::new(name, seed)));
                &current.as_ref().expect("just set").1
            }
        };
        for position in span.start.get()..=span.end.get() {
            if hasher.hash_position(position) < threshold {
                kept.push(GenomePosition {
                    contig: span.contig,
                    position: Position(position),
                });
            }
        }
    }
    Ok(kept)
}

// ---------------------------------------------------------------------
// What a run was asked for, and what it produced
// ---------------------------------------------------------------------

/// The reference this run selected against, as one comparable value.
///
/// The whole-reference MD5 [`ReferenceInfo`] already computes — every contig's uppercased
/// bases concatenated in file order — so it pins the contig *order* as well as the bases,
/// which matters here because [`hash_position`] is keyed by a contig's name and two
/// references holding the same contigs in a different order would select different
/// positions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ReferenceDigest(pub [u8; 16]);

impl ReferenceDigest {
    /// The digest of a reference that was read from a FASTA.
    ///
    /// # Errors
    ///
    /// [`SelectionError::ReferenceNotDigested`] when the reference was read from a `.fai`
    /// alone. A `.fai` describes a genome's geometry and holds no bases, so there is
    /// nothing to digest — and a selection whose terms say nothing about the bases it
    /// selected over cannot refuse a mismatched run.
    pub fn of(reference: &ReferenceInfo) -> Result<Self, SelectionError> {
        reference
            .md5
            .map(Self)
            .ok_or(SelectionError::ReferenceNotDigested)
    }
}

/// The analysed region set, as one comparable value.
///
/// **The likeliest of the seven to differ by accident**, because a BED feels like a runtime
/// convenience and is not: two samples walked under different BEDs selected different loci
/// from different denominators, and every other value in [`SelectionTerms`] agrees.
///
/// Digests the spans as coordinates rather than the file's bytes, so two BEDs that name the
/// same territory with different whitespace, ordering or track lines compare equal — the
/// question is which positions were in play, not which file said so.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RegionSetDigest(pub [u8; 16]);

impl RegionSetDigest {
    pub fn of(regions: &SelectableRegions) -> Self {
        let mut hasher = Md5::new();
        for span in regions.spans() {
            hasher.update(span.contig.get().to_le_bytes());
            hasher.update(span.start.get().to_le_bytes());
            hasher.update(span.end.get().to_le_bytes());
        }
        Self(hasher.finalize().into())
    }
}

/// What the catalog file on disk was built at.
///
/// **A different weighting is a different set of tracts, not a subset of one**, so this is
/// compared for equality rather than for permissiveness — and `tool_version` is in it
/// because a change in the detector invalidates the file even when every setting matches.
#[derive(Debug, Clone)]
pub struct CatalogBuildSettings {
    pub criteria: StrRepeatCriteria,
    pub scan: ScanParams,
    pub tool_version: String,
}

impl CatalogBuildSettings {
    pub fn of(catalog: &RepeatCatalog) -> Self {
        let header = catalog.header();
        Self {
            criteria: header.built_under.clone(),
            scan: header.scan,
            tool_version: header.tool_version.clone(),
        }
    }
}

impl PartialEq for CatalogBuildSettings {
    /// **Destructured without `..` on purpose**: a field added to this struct stops this
    /// function compiling instead of quietly dropping out of the comparison, and a value
    /// that drops out of the comparison turns a refusal to pool two runs into a silent
    /// pooling of evidence recorded under different settings.
    fn eq(&self, other: &Self) -> bool {
        let Self {
            criteria,
            scan,
            tool_version,
        } = self;
        let Self {
            criteria: their_criteria,
            scan: their_scan,
            tool_version: their_tool_version,
        } = other;
        same_criteria(criteria, their_criteria)
            && scan == their_scan
            && tool_version == their_tool_version
    }
}

/// Two sets of STR criteria, compared as whole values rather than field by field — a
/// difference anywhere means the two runs asked the catalog for different loci.
///
/// **The purity floor is compared by bit pattern.** `StrRepeatCriteria` carries an `f32`
/// purity floor, and the `==` that `derive(PartialEq)` generates answers `false` for two
/// `NaN` floors that came from the same configuration file — which is the one answer this
/// check must not give, since it would refuse a run that is in fact identical.
fn same_criteria(left: &StrRepeatCriteria, right: &StrRepeatCriteria) -> bool {
    if left.classification.min_purity.to_bits() != right.classification.min_purity.to_bits() {
        return false;
    }
    // Everything else compares by value. The float is neutralised in both copies first, so
    // this comparison cannot reach it — and adding a field to `StrRepeatCriteria` keeps
    // being covered here without anyone remembering to come back.
    let (mut left, mut right) = (left.clone(), right.clone());
    left.classification.min_purity = 0.0;
    right.classification.min_purity = 0.0;
    left == right
}

/// Everything that decides which loci a run keeps — the seven values that cross the machine
/// boundary.
///
/// The samples are walked separately and never see each other's choices, so this is what
/// lets the fit **refuse** to pool two runs rather than average them. There is no tolerance
/// and no partial match: a set difference is meaningless rather than noisy.
///
/// **`PartialEq` and not `Eq`, and that is forced rather than chosen** — see
/// [`same_criteria`].
#[derive(Debug, Clone)]
pub struct SelectionTerms {
    /// The run's selection seed. A different seed selects a disjoint set, and both sets
    /// look well-formed.
    pub seed: u64,
    /// Which reference.
    pub reference: ReferenceDigest,
    /// Which analysed regions.
    pub analysed_regions: RegionSetDigest,
    /// What the catalog file was built at.
    pub catalog_built_under: CatalogBuildSettings,
    /// What *this run* asked the catalog for — a reader chooses its floors freely within
    /// what the file was built at, so this is a second value and not the one above.
    pub ssr_criteria: StrRepeatCriteria,
    /// The generic target position count.
    pub generic_target: u64,
    /// The STR per-stratum cap.
    pub ssr_cap: usize,
}

impl PartialEq for SelectionTerms {
    /// **Destructured without `..` on purpose** — see [`CatalogBuildSettings::eq`]. These are
    /// the values that let the fit refuse two runs walked on different machines, so a field
    /// that silently left the comparison would pool incomparable evidence and never panic.
    fn eq(&self, other: &Self) -> bool {
        let Self {
            seed,
            reference,
            analysed_regions,
            catalog_built_under,
            ssr_criteria,
            generic_target,
            ssr_cap,
        } = self;
        seed == &other.seed
            && reference == &other.reference
            && analysed_regions == &other.analysed_regions
            && catalog_built_under == &other.catalog_built_under
            && same_criteria(ssr_criteria, &other.ssr_criteria)
            && generic_target == &other.generic_target
            && ssr_cap == &other.ssr_cap
    }
}

impl SelectionTerms {
    /// Which field two sets of selection terms first disagree on, in the order a reader
    /// would check them — `None` when they agree.
    ///
    /// **The fit reports this rather than "the terms differ"**, because every value here
    /// fails the same way and only the name says what to fix.
    ///
    /// **Destructured without `..` on purpose** — see [`SelectionTerms::eq`]. A field added
    /// to this struct stops this function compiling rather than silently going unchecked.
    pub fn first_disagreement(&self, other: &Self) -> Option<&'static str> {
        let Self {
            seed,
            reference,
            analysed_regions,
            catalog_built_under,
            ssr_criteria,
            generic_target,
            ssr_cap,
        } = self;
        if seed != &other.seed {
            return Some("selection seed");
        }
        if reference != &other.reference {
            return Some("reference digest");
        }
        if analysed_regions != &other.analysed_regions {
            return Some("analysed region set");
        }
        if catalog_built_under != &other.catalog_built_under {
            return Some("repeat catalog build settings");
        }
        if !same_criteria(ssr_criteria, &other.ssr_criteria) {
            return Some("STR routing criteria");
        }
        if generic_target != &other.generic_target {
            return Some("generic target position count");
        }
        if ssr_cap != &other.ssr_cap {
            return Some("STR per-stratum cap");
        }
        None
    }
}

/// The loci every sample keeps raw evidence at, for one run.
///
/// **Reproducible, not transported**: two samples that build this from the same
/// [`SelectionTerms`] get equal values, which is what lets samples be walked on different
/// machines. Nothing here is written to a sample's records — [`SelectionTerms`] is, and
/// [`CensusLociDigest`] witnesses that the rule really did produce the same list on both.
pub struct CensusLoci {
    generic: Vec<GenomePosition>,
    ssr: StratumSample,
    ssr_stratum_counts: StratumCounts,
}

impl CensusLoci {
    /// Assemble from parts.
    ///
    /// **[`select_kept_loci`] is how a run gets one.** This exists for the consumer that
    /// already holds the three pieces — a fit reading a stored selection back, and the tests
    /// that need a selection without a catalog behind it.
    pub fn from_parts(
        generic: Vec<GenomePosition>,
        ssr: StratumSample,
        ssr_stratum_counts: StratumCounts,
    ) -> Self {
        Self {
            generic,
            ssr,
            ssr_stratum_counts,
        }
    }

    /// The generic set, ascending in genome order.
    ///
    /// **The order is part of the contract**, not an artifact of how it was built: the fit
    /// indexes a sample's records by position in this slice, and the records carry no
    /// coordinates of their own.
    pub fn generic(&self) -> &[GenomePosition] {
        &self.generic
    }

    /// The STR set, keyed by *(motif period, reference repeat count)*.
    pub fn ssr(&self) -> &StratumSample {
        &self.ssr
    }

    /// How many loci each stratum **holds**, against the [`ssr`](Self::ssr) it kept.
    ///
    /// Anything pooled across strata is biased without this and silently so, which is why
    /// it travels beside the sample rather than being recoverable from it.
    pub fn ssr_stratum_counts(&self) -> &StratumCounts {
        &self.ssr_stratum_counts
    }
}

/// A witness that two samples really did keep the same loci.
///
/// The seven values of [`SelectionTerms`] say what a run was *asked*; this says what it
/// *produced*. They fail differently: the seven all agree when a hash function or a
/// threshold's arithmetic has changed underneath them, and this does not.
///
/// **It must be fed as records are written**, one call per kept locus in index order, so it
/// digests the array that exists rather than the list that should have been built. A digest
/// derived by re-running the selection proves only that the selection is deterministic.
///
/// **Blocked per megabase so a mismatch names where it happened**, and each block carries
/// the coordinates it covers rather than only a number — 800 blocks on tomato at a
/// two-million-position budget, which is nothing beside the records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CensusLociDigest {
    whole: [u8; 16],
    per_block: Vec<BlockDigest>,
}

/// One megabase of the analysed regions, and the digest of the loci kept inside it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BlockDigest {
    pub contig: ContigId,
    /// Which megabase of that contig — the position divided by 1,000,000.
    pub megabase: u32,
    pub digest: u64,
}

/// Fills a [`CensusLociDigest`] as the [`CensusWriter`](super::census::CensusWriter) walks
/// the kept loci.
///
/// Separate from the finished digest so that the finished one is immutable and comparable,
/// and so `observe`'s index check has somewhere to keep its expectation.
#[derive(Default)]
pub struct CensusLociDigester {
    whole: Md5,
    blocks: Vec<BlockDigest>,
    open_block: Option<(ContigId, u32, Xxh3)>,
    next_index: usize,
}

/// A megabase, in bases. The block width the per-block digests are cut at.
const DIGEST_BLOCK_BP: u64 = 1_000_000;

impl CensusLociDigester {
    pub fn new() -> Self {
        Self::default()
    }

    /// Absorb the `index`-th kept locus.
    ///
    /// # Panics
    ///
    /// When `index` is not the one expected next. **A writer that skips or reorders must
    /// fail loudly**: the digest it would otherwise produce is plausible, matches nothing,
    /// and names no reason.
    pub fn observe(&mut self, index: usize, locus: GenomePosition) {
        assert_eq!(
            index, self.next_index,
            "kept loci must be digested in index order; a skipped or reordered locus \
             produces a plausible digest that witnesses nothing"
        );
        self.next_index += 1;

        self.whole.update(locus.contig.get().to_le_bytes());
        self.whole.update(locus.position.get().to_le_bytes());

        let megabase = u32::try_from(locus.position.get() / DIGEST_BLOCK_BP)
            .expect("a contig shorter than four petabases");
        let reopen = match &self.open_block {
            Some((contig, open, _)) => *contig != locus.contig || *open != megabase,
            None => true,
        };
        if reopen {
            self.close_block();
            self.open_block = Some((locus.contig, megabase, Xxh3::new()));
        }
        let (_, _, hasher) = self.open_block.as_mut().expect("just opened");
        locus.position.get().hash(hasher);
    }

    fn close_block(&mut self) {
        if let Some((contig, megabase, hasher)) = self.open_block.take() {
            self.blocks.push(BlockDigest {
                contig,
                megabase,
                digest: hasher.finish(),
            });
        }
    }

    /// How many loci have been absorbed — the count the caller checks against
    /// [`CensusLoci::generic`]'s length before trusting the digest.
    pub fn observed(&self) -> usize {
        self.next_index
    }

    pub fn finish(mut self) -> CensusLociDigest {
        self.close_block();
        CensusLociDigest {
            whole: self.whole.finalize().into(),
            per_block: self.blocks,
        }
    }
}

impl CensusLociDigest {
    /// The first megabase two digests disagree on, as `(contig, megabase)` — `None` when
    /// they agree.
    ///
    /// **This is the whole reason for blocking.** Two samples whose whole digests differ
    /// have kept different loci somewhere, and "somewhere" in two million positions is not
    /// a thing anyone can act on.
    pub fn first_disagreement(&self, other: &Self) -> Option<(ContigId, u32)> {
        for (mine, theirs) in self.per_block.iter().zip(&other.per_block) {
            if mine != theirs {
                return Some((mine.contig, mine.megabase));
            }
        }
        // One ran out first: the first block the shorter one does not have.
        let shared = self.per_block.len().min(other.per_block.len());
        let longer = if self.per_block.len() > other.per_block.len() {
            &self.per_block
        } else {
            &other.per_block
        };
        longer
            .get(shared)
            .map(|block| (block.contig, block.megabase))
    }

    pub fn blocks(&self) -> &[BlockDigest] {
        &self.per_block
    }
}

// ---------------------------------------------------------------------
// The selection itself
// ---------------------------------------------------------------------

/// Derive the kept loci from the run's inputs alone — **no read is opened**.
///
/// `analysed` is the run's own territory (every contig whole with no `--regions` BED);
/// `unambiguous` is [`UnambiguousRuns::into_selectable`]'s output for the same reference.
/// The two are intersected here rather than by the caller, because the generic threshold's
/// denominator has to be the same object the sweep walks and separating them is the failure
/// this module is shaped to prevent.
///
/// **The generic domain excludes repeat tracts as well as ambiguous bases.** A tract's
/// variability would distort every substitution statistic computed from the set, and the STR
/// path is where tracts are kept — so the domain is the analysed regions, masked, then cut
/// by the catalog and reduced to the [`Generic`](RegionKind::Generic) pieces. The threshold
/// is computed from *that* length.
///
/// **The STR scope is the analysed regions unmasked**, because a tract in the catalog is a
/// stretch of real sequence by construction, and cutting the scope at every `N` run would
/// split scope spans without changing which tracts they hold.
///
/// # Errors
///
/// Passes through the catalog's refusal when the run asks for a copy floor or a flank the
/// file was not built at — the catalog refuses rather than serving a short list, and a short
/// list here is a stratum that looks depleted rather than unasked for.
pub fn select_kept_loci(
    terms: &SelectionTerms,
    catalog: &RepeatCatalog,
    analysed: &SelectableRegions,
    unambiguous: &SelectableRegions,
) -> Result<CensusLoci, SelectionError> {
    let names: Vec<String> = catalog.contigs().iter().map(|c| c.name.clone()).collect();

    let masked = analysed.intersect(unambiguous);
    let mut generic_spans = Vec::new();
    for segment in catalog
        .genome_segments(&terms.ssr_criteria, ReadScope::Regions(masked.spans()))
        .map_err(SelectionError::Catalog)?
    {
        let segment = segment.map_err(SelectionError::Catalog)?;
        if segment.kind == RegionKind::Generic {
            generic_spans.push(segment.region);
        }
    }
    let generic_domain = SelectableRegions::new(generic_spans)?;

    let threshold = threshold_for(terms.generic_target, generic_domain.total_length());
    let generic = select_generic_positions(&generic_domain, &names, terms.seed, threshold)?;

    let (ssr_stratum_counts, ssr) = catalog
        .sample_loci_per_stratum(
            &terms.ssr_criteria,
            ReadScope::Regions(analysed.spans()),
            terms.ssr_cap,
            terms.seed,
        )
        .map_err(SelectionError::Catalog)?;

    Ok(CensusLoci {
        generic,
        ssr,
        ssr_stratum_counts,
    })
}

// ---------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------

/// Why a selection could not be built. **Every variant is fatal**: a selection that is
/// wrong is not a degraded estimate but a meaningless one, since every sample would have
/// to be wrong in the identical way for the fit to notice nothing.
///
/// **Not `PartialEq`**, because [`Catalog`](Self::Catalog) carries the catalog's own
/// refusal and that type is not comparable. A test asks `matches!` rather than `assert_eq!`
/// — which is the right question anyway: what a caller acts on is which variant, not which
/// coordinates were in the message.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum SelectionError {
    #[error("selectable region {region:?} ends before it starts")]
    InvertedRegion { region: GenomeRegion },
    #[error("selectable regions {first:?} and {second:?} overlap; a position would be kept twice")]
    OverlappingRegions {
        first: GenomeRegion,
        second: GenomeRegion,
    },
    #[error("selectable regions name contig {contig:?}, which the reference does not hold")]
    UnknownContig { contig: ContigId },
    /// The reference was read from a `.fai` alone, so there are no bases to digest and the
    /// selection terms could not say which reference this selection was made over.
    #[error("the reference was read without its bases, so it has no content digest")]
    ReferenceNotDigested,
    #[error("the repeat catalog cannot serve this selection")]
    Catalog(#[source] RepeatCatalogError),
}

#[cfg(test)]
mod tests {
    use super::*;

    fn names() -> Vec<String> {
        vec!["chr1".to_string(), "chr2".to_string()]
    }

    fn region(contig: u32, start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(contig),
            start: Position(start),
            end: Position(end),
        }
    }

    #[test]
    fn hoisting_the_contig_name_changes_no_hash() {
        let hasher = ContigHasher::new("chr1", 42);
        for position in [1_u64, 2, 1_000, u32::MAX as u64, u64::MAX] {
            assert_eq!(
                hasher.hash_position(position),
                hash_position("chr1", position, 42),
                "the hoisted hasher must agree with the one-off form at {position}"
            );
        }
    }

    #[test]
    fn the_threshold_is_computed_in_128_bits() {
        // 2^64 · 1 / 4 — the arithmetic that overflows a u64 if done there.
        assert_eq!(threshold_for(1, Bp(4)), 1_u64 << 62);
        assert_eq!(threshold_for(1, Bp(2)), 1_u64 << 63);
        // A target at or above the length keeps everything rather than being an error.
        assert_eq!(threshold_for(10, Bp(10)), u64::MAX);
        assert_eq!(threshold_for(11, Bp(10)), u64::MAX);
        assert_eq!(threshold_for(0, Bp(10)), 0);
        // An empty domain has nothing to select from and must not divide by zero.
        assert_eq!(threshold_for(5, Bp(0)), u64::MAX);
    }

    #[test]
    fn the_realised_count_lands_near_the_target() {
        let regions = SelectableRegions::new(vec![region(0, 1, 1_000_000)]).unwrap();
        let target = 10_000;
        let threshold = threshold_for(target, regions.total_length());
        let kept = select_generic_positions(&regions, &names(), 7, threshold).unwrap();
        let error = (kept.len() as i64 - target as i64).abs();
        // Binomial around the target: ±4√target is a 4-sigma band, so this fails on a
        // broken hash and not on an unlucky seed.
        assert!(
            error < 4 * (target as f64).sqrt() as i64,
            "kept {} against a target of {target}",
            kept.len()
        );
    }

    #[test]
    fn a_sharded_sweep_keeps_the_same_positions() {
        let whole =
            SelectableRegions::new(vec![region(0, 1, 200_000), region(1, 1, 200_000)]).unwrap();
        let threshold = threshold_for(2_000, whole.total_length());
        let from_whole = select_generic_positions(&whole, &names(), 11, threshold).unwrap();

        // The same domain cut into shards that do not respect any boundary the rule knows
        // about. The threshold is the run's, not the shard's — a shard that recomputed it
        // from its own length would keep the right *fraction* of a different denominator.
        let sharded = SelectableRegions::new(vec![
            region(0, 1, 37),
            region(0, 38, 149_999),
            region(0, 150_000, 200_000),
            region(1, 1, 100_000),
            region(1, 100_001, 200_000),
        ])
        .unwrap();
        let from_shards = select_generic_positions(&sharded, &names(), 11, threshold).unwrap();
        assert_eq!(from_whole, from_shards);
    }

    #[test]
    fn kept_positions_are_sorted_and_inside_the_regions() {
        let regions =
            SelectableRegions::new(vec![region(1, 5_000, 9_999), region(0, 100, 50_000)]).unwrap();
        let threshold = threshold_for(500, regions.total_length());
        let kept = select_generic_positions(&regions, &names(), 3, threshold).unwrap();
        assert!(
            kept.len() > 100,
            "too few kept to be a test: {}",
            kept.len()
        );
        for pair in kept.windows(2) {
            assert!(
                (pair[0].contig.get(), pair[0].position.get())
                    < (pair[1].contig.get(), pair[1].position.get()),
                "kept positions must ascend and never repeat: {pair:?}"
            );
        }
        for position in &kept {
            assert!(
                regions
                    .spans()
                    .iter()
                    .any(|span| span.contig == position.contig
                        && span.start.get() <= position.position.get()
                        && position.position.get() <= span.end.get()),
                "{position:?} lies outside the selectable regions"
            );
        }
    }

    #[test]
    fn a_different_seed_selects_a_different_set() {
        let regions = SelectableRegions::new(vec![region(0, 1, 100_000)]).unwrap();
        let threshold = threshold_for(1_000, regions.total_length());
        let one = select_generic_positions(&regions, &names(), 1, threshold).unwrap();
        let two = select_generic_positions(&regions, &names(), 2, threshold).unwrap();
        assert_ne!(one, two);
        let shared = one.iter().filter(|p| two.contains(p)).count();
        // Two independent draws of one position in a hundred overlap in about one
        // position in ten thousand — a handful out of a thousand, never most of them.
        assert!(
            shared * 20 < one.len(),
            "{shared} of {} positions shared between two seeds",
            one.len()
        );
    }

    #[test]
    fn the_selection_is_narrowed_by_the_region_set_and_so_is_the_denominator() {
        // The same target over a tenth of the genome keeps the target, not a tenth of it:
        // the threshold reads the region set's own length.
        let whole = SelectableRegions::new(vec![region(0, 1, 1_000_000)]).unwrap();
        let tenth = SelectableRegions::new(vec![region(0, 1, 100_000)]).unwrap();
        let from_tenth = select_generic_positions(
            &tenth,
            &names(),
            5,
            threshold_for(5_000, tenth.total_length()),
        )
        .unwrap();
        assert_eq!(whole.total_length(), Bp(1_000_000));
        assert_eq!(tenth.total_length(), Bp(100_000));
        let error = (from_tenth.len() as i64 - 5_000).abs();
        assert!(error < 4 * 71, "kept {} against 5,000", from_tenth.len());
    }

    #[test]
    fn overlapping_regions_are_refused() {
        let err = SelectableRegions::new(vec![region(0, 1, 100), region(0, 50, 200)]).unwrap_err();
        assert!(matches!(err, SelectionError::OverlappingRegions { .. }));
        // Abutting is not overlapping.
        assert!(SelectableRegions::new(vec![region(0, 1, 100), region(0, 101, 200)]).is_ok());
        // The same coordinates on two contigs are two different places.
        assert!(SelectableRegions::new(vec![region(0, 1, 100), region(1, 1, 100)]).is_ok());
    }

    #[test]
    fn the_unambiguous_mask_keeps_the_runs_and_drops_the_gaps() {
        let mut mask = UnambiguousRuns::default();
        mask.contig_started("chr1", 0);
        // Fed in three batches, so a run that spans a batch boundary has to survive one.
        mask.bases(b"ACGTNN");
        mask.bases(b"NACG");
        mask.bases(b"TnACGT");
        mask.contig_finished(&ContigInfo {
            name: "chr1".to_string(),
            length: 16,
            offset: 0,
            line_bases: 16,
            line_width: 17,
            md5: None,
        });
        assert_eq!(mask.ambiguous_bases(), 4);
        let regions = mask.into_selectable().unwrap();
        assert_eq!(
            regions.spans(),
            &[region(0, 1, 4), region(0, 8, 11), region(0, 13, 16)]
        );
        // The threshold's denominator is the mask's length, not the contig's.
        assert_eq!(regions.total_length(), Bp(12));
    }

    #[test]
    fn a_span_naming_an_absent_contig_stops_the_run() {
        let regions = SelectableRegions::new(vec![region(9, 1, 100)]).unwrap();
        let err = select_generic_positions(&regions, &names(), 1, u64::MAX).unwrap_err();
        assert!(matches!(
            err,
            SelectionError::UnknownContig {
                contig: ContigId(9)
            }
        ));
    }

    // ---- narrowing the domain ------------------------------------------------

    #[test]
    fn intersecting_narrows_the_spans_and_the_length_together() {
        // Two analysed stretches; the mask cuts a hole in the first and clips the second.
        let analysed =
            SelectableRegions::new(vec![region(0, 1, 100), region(0, 200, 300)]).unwrap();
        let unambiguous =
            SelectableRegions::new(vec![region(0, 50, 249), region(0, 260, 1000)]).unwrap();

        let masked = analysed.intersect(&unambiguous);

        assert_eq!(
            masked.spans(),
            [region(0, 50, 100), region(0, 200, 249), region(0, 260, 300)]
        );
        // 51 + 50 + 41, and the length is the sum of what is there rather than a separate
        // claim about it — the failure the type exists to prevent.
        assert_eq!(masked.total_length(), Bp(142));
        assert_eq!(
            masked.total_length().get(),
            masked
                .spans()
                .iter()
                .map(|s| s.end.get() - s.start.get() + 1)
                .sum::<u64>()
        );
    }

    #[test]
    fn intersecting_keeps_contigs_apart() {
        let analysed = SelectableRegions::new(vec![region(0, 1, 100), region(1, 1, 100)]).unwrap();
        let unambiguous = SelectableRegions::new(vec![region(1, 40, 60)]).unwrap();
        let masked = analysed.intersect(&unambiguous);
        assert_eq!(masked.spans(), [region(1, 40, 60)]);
    }

    // ---- the selection terms -------------------------------------------------

    fn terms() -> SelectionTerms {
        SelectionTerms {
            seed: 42,
            reference: ReferenceDigest([7; 16]),
            analysed_regions: RegionSetDigest([9; 16]),
            catalog_built_under: CatalogBuildSettings {
                criteria: StrRepeatCriteria::default(),
                scan: ScanParams::default(),
                tool_version: "0.1.0".to_string(),
            },
            ssr_criteria: StrRepeatCriteria::default(),
            generic_target: 2_000_000,
            ssr_cap: 1_000,
        }
    }

    #[test]
    fn a_set_of_terms_equals_itself_even_with_a_not_a_number_purity_floor() {
        // The one case `derive(PartialEq)` gets wrong: a run refusing itself because a
        // configuration file put a NaN in the purity floor of both copies.
        let mut left = terms();
        left.ssr_criteria.classification.min_purity = f32::NAN;
        let right = left.clone();
        assert_eq!(left, right);
        assert_eq!(left.first_disagreement(&right), None);
    }

    #[test]
    fn a_different_purity_floor_is_a_different_selection() {
        let left = terms();
        let mut right = left.clone();
        right.ssr_criteria.classification.min_purity += 0.01;
        assert_ne!(left, right);
        assert_eq!(
            left.first_disagreement(&right),
            Some("STR routing criteria")
        );
    }

    #[test]
    fn each_field_is_named_when_it_is_the_one_that_differs() {
        let left = terms();

        let mut seed = left.clone();
        seed.seed += 1;
        assert_eq!(left.first_disagreement(&seed), Some("selection seed"));

        let mut reference = left.clone();
        reference.reference = ReferenceDigest([0; 16]);
        assert_eq!(
            left.first_disagreement(&reference),
            Some("reference digest")
        );

        let mut regions = left.clone();
        regions.analysed_regions = RegionSetDigest([0; 16]);
        assert_eq!(
            left.first_disagreement(&regions),
            Some("analysed region set")
        );

        let mut built = left.clone();
        built.catalog_built_under.tool_version = "0.2.0".to_string();
        assert_eq!(
            left.first_disagreement(&built),
            Some("repeat catalog build settings")
        );

        let mut target = left.clone();
        target.generic_target += 1;
        assert_eq!(
            left.first_disagreement(&target),
            Some("generic target position count")
        );

        let mut cap = left.clone();
        cap.ssr_cap += 1;
        assert_eq!(left.first_disagreement(&cap), Some("STR per-stratum cap"));
    }

    #[test]
    fn two_region_sets_naming_the_same_territory_digest_alike() {
        // Same spans, offered in a different order: `SelectableRegions::new` sorts, so the
        // digest is over the territory rather than over how it was written down.
        let one = SelectableRegions::new(vec![region(0, 1, 100), region(1, 5, 50)]).unwrap();
        let other = SelectableRegions::new(vec![region(1, 5, 50), region(0, 1, 100)]).unwrap();
        assert_eq!(RegionSetDigest::of(&one), RegionSetDigest::of(&other));

        let shifted = SelectableRegions::new(vec![region(0, 1, 101), region(1, 5, 50)]).unwrap();
        assert_ne!(RegionSetDigest::of(&one), RegionSetDigest::of(&shifted));
    }

    // ---- the digest that checks the answer -----------------------------------

    fn locus(contig: u32, position: u64) -> GenomePosition {
        GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        }
    }

    fn digest_of(loci: &[GenomePosition]) -> CensusLociDigest {
        let mut digester = CensusLociDigester::new();
        for (index, locus) in loci.iter().enumerate() {
            digester.observe(index, *locus);
        }
        assert_eq!(digester.observed(), loci.len());
        digester.finish()
    }

    #[test]
    fn the_same_loci_digest_alike_and_one_moved_locus_does_not() {
        let loci = vec![locus(0, 10), locus(0, 1_500_000), locus(1, 40)];
        assert_eq!(digest_of(&loci), digest_of(&loci));

        let mut moved = loci.clone();
        moved[1] = locus(0, 1_500_001);
        assert_ne!(digest_of(&loci), digest_of(&moved));
    }

    #[test]
    fn a_disagreement_names_the_megabase_it_happened_in() {
        let loci = vec![locus(0, 10), locus(0, 1_500_000), locus(1, 40)];
        let mut moved = loci.clone();
        moved[1] = locus(0, 1_500_001);

        // The first block — contig 0's megabase 0 — is untouched, so the answer must be the
        // second, not merely "they differ".
        assert_eq!(
            digest_of(&loci).first_disagreement(&digest_of(&moved)),
            Some((ContigId(0), 1))
        );
        assert_eq!(digest_of(&loci).first_disagreement(&digest_of(&loci)), None);
    }

    #[test]
    fn a_selection_that_stops_early_is_caught_by_the_block_it_never_reached() {
        let full = vec![locus(0, 10), locus(0, 1_500_000), locus(1, 40)];
        let short = &full[..2];
        assert_eq!(
            digest_of(&full).first_disagreement(&digest_of(short)),
            Some((ContigId(1), 0))
        );
    }

    #[test]
    fn one_block_per_megabase_occupied() {
        let loci = vec![
            locus(0, 10),
            locus(0, 999_999),
            locus(0, 1_000_000),
            locus(1, 5),
        ];
        let blocks = digest_of(&loci);
        let named: Vec<_> = blocks
            .blocks()
            .iter()
            .map(|b| (b.contig, b.megabase))
            .collect();
        assert_eq!(
            named,
            [(ContigId(0), 0), (ContigId(0), 1), (ContigId(1), 0)]
        );
    }

    #[test]
    #[should_panic(expected = "index order")]
    fn digesting_out_of_order_fails_loudly() {
        let mut digester = CensusLociDigester::new();
        digester.observe(0, locus(0, 10));
        digester.observe(2, locus(0, 20));
    }
}
