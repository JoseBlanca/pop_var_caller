//! ng step 3 — what the sequence at a place in the reference *is*: the vocabulary of
//! typed regions, the policy that decides between them, and one implementation that
//! answers the question straight from a contig's bases.
//! Design: `doc/devel/ng/spec/typed_regions.md` (spec) and
//! `doc/devel/ng/arch/typed_regions.md` (types & interfaces).
//!
//! **A run does not answer this question here.** Every consumer reads the reference's
//! repeat catalog — the tandem repeats found once per reference and written beside the
//! FASTA — and derives its regions from that file
//! ([`crate::ng::repeat_catalog`], spec `repeat_catalog.md`). What this module still owns
//! is everything after detection: the region vocabulary ([`TypedRegion`], [`RegionKind`]),
//! the policy ([`TypedRegionConfig`], [`segment_criteria`]), the satellite cap, the
//! bundling, the generic fill and the tally ([`TypedRegionCounts`]) — all of which the
//! catalog's reader calls rather than copies.
//!
//! [`partition_resident`] is the one thing here that still opens with a scan, and it
//! exists to be the catalog's yardstick: it derives the same segmentation from the bases,
//! so the differential can say the file gives the same answer rather than a similar one.
//!
//! **A folder, not a file, and not because of a bake-off** (there is none —
//! spec §6). The classification port is a second concern with its own dense test
//! suite, so it gets its own module.
//!
//! ## Production is frozen; ng owns its copies
//!
//! Step 3 needs an STR classification policy that is 1-based/`u64`, driven by
//! `RepeatInterval`s, all-knobs, and that hands bundle members back instead of dropping
//! them. `ssr::catalog::postprocess::build_loci` is none of those things, and **reshaping
//! it in place is not on the table** (spec Revision 2026-07-16, owner): production stays
//! exactly as it is, so that it remains an *independent yardstick* for the experiments ng
//! exists to run.
//!
//! So [`segment_criteria`] is a **port**: the logic transcribed unchanged, the shape
//! ng's. What sharing one function used to guarantee for free, a test now pins
//! — see [`segment_criteria`]'s differential against production (spec §8.0).

pub mod segment_criteria;

use std::path::Path;

use crate::ng::tandem_repeat::{RepeatInterval, ScanParams, find_tandem_repeats};
use crate::ng::types::{Bp, ContigId, GenomeRegion, Position};
use crate::regions::{BedError, ContigBounds, RegionSet};
use segment_criteria::{RejectionCounts, SsrSegment, SsrSegmentCriteria};

// ---------------------------------------------------------------------
// What to walk
// ---------------------------------------------------------------------

/// The set of genome regions to walk — sorted, non-overlapping, coalesced,
/// clamped, in genomic order.
///
/// **Wraps production's `RegionSet` read-only; reimplements nothing.** That type
/// already parses BED, coalesces overlapping and adjacent spans, clamps to contig
/// lengths, resolves names against the contig table, and drops zero-length
/// contigs — and it is the same code the production caller's `--regions` runs, so
/// ng and production agree on what a BED *means* by construction rather than by
/// coincidence. `src/regions.rs` is not edited (spec Revision): this wrapper adds
/// ng's width and ng's names, and nothing else.
///
/// **A user BED is not a special case**, which `regions.rs` settled first: *"'Whole
/// genome' is not a special case — it is the region set whose every region covers
/// an entire contig."* [`Self::whole_contigs`] is the default, not a bypass.
///
/// ## The conversion this owns — and it is smaller than the spec expected
///
/// Spec §4 called this "the one conversion seam", *"`GenomeRegions` widening (and
/// **rebasing**) `RegionSet`'s `u32`"*. **There is no rebasing.** `regions::Region`
/// is already **1-based inclusive** — its own doc says so, and its invariant is
/// `1 <= start <= end`. So production and ng already agree on the base, and the
/// only conversion here is `u32` → `u64`, which is lossless and infallible.
///
/// That is worth stating rather than quietly enjoying: the spec anticipated an
/// off-by-one seam at ng's busiest boundary and there is none, because the
/// production author had already made the same call for the same reason.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenomeRegions {
    inner: RegionSet,
}

impl GenomeRegions {
    /// One full-length span per contig — **the default** (spec §2.5).
    ///
    /// Zero-length contigs contribute no span, so they never reach a consumer
    /// (`RegionSet`'s rule, and the reason spec §2.3 can say "zero-length contigs
    /// are never asked about" without a guard of its own).
    pub fn whole_contigs(contigs: &[ContigBounds]) -> Self {
        Self {
            inner: RegionSet::whole_contigs(contigs),
        }
    }

    /// Parse a BED and resolve it against the contig table.
    ///
    /// Every failure mode — a short line, non-numeric coordinates, `end <= start`,
    /// an unknown contig name, a span past a contig's end — is `RegionSet`'s to
    /// reject, **up front** (spec §8.2). So a consumer holding one of these has
    /// nothing left to validate: the spans name contigs the reference has and stay
    /// inside them.
    pub fn from_bed_path(bed: &Path, contigs: &[ContigBounds]) -> Result<Self, BedError> {
        Ok(Self {
            inner: RegionSet::from_bed_path(bed, contigs)?,
        })
    }

    /// The regions, in genomic order, as ng's [`GenomeRegion`].
    ///
    /// The `u32` → `u64` widening lives here and only here (above).
    pub fn iter(&self) -> impl Iterator<Item = GenomeRegion> + '_ {
        self.inner.iter().map(|r| GenomeRegion {
            contig: ContigId(r.chrom_id),
            start: Position(u64::from(r.start)),
            end: Position(u64::from(r.end)),
        })
    }

    /// How many regions there are.
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    /// Whether there is nothing to ask about.
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

// ---------------------------------------------------------------------
// What comes out
// ---------------------------------------------------------------------

/// A genome region plus **what the sequence there is**.
///
/// `region` is a field, not a per-variant repeat: every typed region has one,
/// *structurally*, and the partition invariant (spec §2.3) reads it off directly.
/// It is also the one place ng's 1-based base is stated for this step (spec §4).
///
/// Rejected: a `region()` accessor over four variants, which makes "every typed
/// region has a region" a convention rather than a fact of the type. The evidence
/// it fails: an earlier draft written against the accessor design stated the
/// invariant in **0-based** phrasing — in the very property the spec calls its
/// spine.
#[derive(Debug, Clone, PartialEq)]
pub struct TypedRegion {
    pub region: GenomeRegion,
    pub kind: RegionKind,
}

/// What the sequence in a region **is** — one of four, and exactly one of them is
/// a *genetic* object (spec §1.1).
///
/// | kind | is | why |
/// |---|---|---|
/// | `SsrSegment` | a **locus** | the only kind the reference alone hands you as a genetic object |
/// | `SsrBundle` | a region | real repeats, none with clean flanks, so no locus can be named |
/// | `Generic` | a region | nothing more specific can be said from the reference alone |
/// | `Satellite` | a region | a tandem array too long to be a microsatellite |
///
/// **This step types regions; it does not decide their fate.** What a consumer
/// then does with each — genotype it, pile it up, mask it, skip it — is a
/// decision downstream (spec §1).
///
/// `Generic` and `Satellite` carry nothing because they *are* just their region.
/// That is not only tidiness: an open generic run costs **two coordinates however
/// many megabases it spans**, which is what makes spec §2.1's emission rule
/// affordable.
#[derive(Debug, Clone, PartialEq)]
pub enum RegionKind {
    /// ng's own [`SsrSegment`] — motif, borders, purity. **Coordinates, no bases**: the
    /// bases are in the reference the caller already has open (see [`SsrSegment`]). No
    /// wrapper: it is 1-based like the rest of ng (spec §4), and [`TypedRegion`]
    /// already carries the region.
    SsrSegment(SsrSegment),
    /// A cluster of repeats none of which has clean flanks (spec §2.4). Carries
    /// the tracts as coordinates — enough to see the structure (each interval has
    /// its period) without this step pre-deciding what it is for. The hull is the
    /// [`TypedRegion`]'s own region.
    ///
    /// `>= 2` members, coordinate-ordered. **This variant is why bundles exist as
    /// a type at all**: production *deletes* these records, so their bases become
    /// a hole nobody accounts for; carrying them lets the decision be taken later,
    /// with the evidence in hand (spec §1, §10).
    SsrBundle { tracts: Box<[RepeatInterval]> },
    /// Nothing more specific can be said from the reference alone — **the
    /// default**, not a leftover. The other three are exceptions carved out of it
    /// (spec §2.2), and a repeat classification turns down for any reason other than
    /// bundling lands back here rather than becoming a hole.
    Generic,
    /// A tandem array longer than `max_str_len` — an array, not a
    /// microsatellite. A **typing** claim, and `max_str_len` is its parameter,
    /// not a constant of nature (spec §2.1, §10).
    Satellite,
}

// ---------------------------------------------------------------------
// Config and counts
// ---------------------------------------------------------------------

/// The policy that decides what the sequence at a place is. Mirrors
/// `ReadFilterConfig`'s shape (`read_filtering.md` §2.4): defaults as named consts,
/// `Default` = what the lab runs, no dormant knobs.
///
/// `Default` is **the catalog's settings, for spec §8's comparability only** — not
/// an endorsement of them (spec §5.2). The catalog is a yardstick, not an
/// authority.
#[derive(Debug, Clone, PartialEq)]
pub struct TypedRegionConfig {
    /// The scanner's scoring weights. **The period range is not here** — it is
    /// [`SsrSegmentCriteria::periods`], the one range the scanner detects and
    /// classification accepts. They were two fields until the scanner learned to emit
    /// only primitive periods (`find_tandem_repeats`); before that the scanner had to
    /// over-scan from period 1 to give the pre-filter an eliminator for homopolymer
    /// aliases, so "detect" and "accept" needed separate ranges. The scanner honours
    /// its range now, so there is one.
    pub scan: ScanParams,
    /// The satellite cap: a tract longer than this is a `Satellite`, not an STR.
    ///
    /// **Also the detection margin of anything that scans in pieces** — the catalog's
    /// builder scans whole contigs and so has none, but `RegionScanner` still does, and
    /// the two must be the same number (spec §2.6): a margin exists to capture whole any
    /// repeat that is not a satellite, so it is exactly the length at which a repeat
    /// becomes one.
    pub max_str_len: Bp,
    /// Admission's rules — all of them (spec §5).
    pub criteria: SsrSegmentCriteria,
}

/// The satellite cap and detection margin: **100 bp** (spec §2.3). This one
/// field is both — a tract longer than the cap is a `Satellite`, not an STR,
/// and a scanner working in pieces fetches core ± this margin. With 150 bp reads a read
/// spans a tract plus an anchor each side only up to ~`read_len − 2·bundle_threshold` ≈ 90 bp,
/// so past ~100 bp the STR route has nothing to offer. A round number at that
/// read-length limit, not a measured one — soft, and the point of the knob is
/// to sweep it (spec §10).
pub const DEFAULT_MAX_STR_LEN: u64 = 100;

impl Default for TypedRegionConfig {
    fn default() -> Self {
        Self {
            scan: ScanParams::default(),
            max_str_len: Bp(DEFAULT_MAX_STR_LEN),
            criteria: SsrSegmentCriteria::default(),
        }
    }
}

/// What a segmentation contained, and what admission turned down getting there
/// ([`partition_resident_in`]).
///
/// **The catalog's [`CatalogRegionCounts`](crate::ng::repeat_catalog::CatalogRegionCounts)
/// is this type minus one counter**, and the two are compared field for field by the
/// differential. Its module doc says which counter is missing and why.
///
/// **"No silent caps"**: a base typed away from the STR path must be accounted
/// for. This is the caller's view of the catalog's measured ~35% STR coverage gap.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct TypedRegionCounts {
    /// Requested spans on this contig — the regions the caller asked about, not the
    /// regions that came back.
    pub spans: u64,
    pub ssr_loci: u64,
    pub ssr_bundles: u64,
    /// **The number spec §10's bundle question needs, and has never had** —
    /// because the answer was previously deleted uncounted (production drops
    /// bundle members without recording them).
    pub ssr_bundle_bp: u64,
    pub generic: u64,
    pub satellites: u64,
    pub satellite_bp: u64,
    /// Repeat coverage that yielded **no locus**, in bp: every base of cleaned repeat
    /// coverage that did not come out as an `SsrSegment` — because it was bundled,
    /// capped as a satellite, or rejected by one of classification's gates.
    ///
    /// In bp, not per repeat, because a per-repeat count answers the wrong
    /// question twice: classification trims every survivor, so a repeat that classifies one
    /// locus and sheds 200 bp contributes nothing to a per-repeat counter (spec
    /// §3.1).
    ///
    /// **Reached by subtraction**: the contig's whole cleaned repeat coverage inside the
    /// request is charged first, and the loci that cancel it are subtracted one by one. The
    /// catalog's reader charges and cancels the same way, which is what keeps the two
    /// numbers one number.
    pub repeat_bp_with_no_locus: u64,
    /// [`Self::repeat_bp_with_no_locus`] **broken out by classification's reason** — because
    /// one total cannot separate a purity rejection from a copy-floor one, and that is
    /// exactly the distinction spec §10's routing question turns on.
    ///
    /// **It does not partition the total**, and [`RejectionCounts`] says why: overlapping
    /// rejected repeats are both charged, and bases with no locus for a reason that is not
    /// a *rejection* (bundled, capped as a satellite, out of scope) are in the total and
    /// not here. A diagnosis of classification's gates, not an account of the genome.
    pub rejected_by_reason: RejectionCounts,
}

// ---------------------------------------------------------------------
// The resident scan — the differential's oracle
// ---------------------------------------------------------------------

/// Cut one whole contig into typed regions, holding it entirely in memory.
///
/// The whole contig asked for, and the tally dropped —
/// [`partition_resident_in`] is the same function with a region subset and a count.
pub fn partition_resident(
    chrom: &str,
    contig: ContigId,
    bases: &[u8],
    config: &TypedRegionConfig,
) -> Vec<TypedRegion> {
    let whole = [GenomeRegion {
        contig,
        start: Position(1),
        end: Position(bases.len().max(1) as u64),
    }];
    partition_resident_in(chrom, contig, bases, config, &whole).0
}

/// The typed regions of one contig that fall inside `wanted`, and what they contained.
///
/// **The five steps** (spec §2.1), and the first three are ng's port of the
/// catalog's implementation (spec §5):
///
/// 1. **Detect** — [`find_tandem_repeats`] → raw, overlapping candidates. No policy.
/// 2. **Clean** — [`segment_criteria::prefilter`]: per-period copy floor, then
///    period-multiple redundancy. **Not optional** (spec §5b).
/// 3. **Admit** — [`segment_criteria::classify`] → the STR loci *and* the tracts it set
///    aside as bundle members.
/// 4. **Cap** — merge the cleaned intervals into coverage runs; a run over
///    `max_str_len` is a satellite. Classified loci inside one are dropped.
/// 5. **Partition** — emit `SsrSegment` at each surviving tract, `Satellite` at each
///    run, `Generic` across everything else, then narrow the result to `wanted`.
///
/// # Why this exists when nothing calls it outside a test
///
/// **It is the catalog's oracle** (`tests/ng_repeat_catalog_differential.rs`). A caller
/// that stops opening the reference and reads the repeat catalog instead is relying on
/// one claim — *the file gives the same segmentation the bases do* — and that claim needs
/// a second implementation to be checked against. This is that implementation: it detects
/// the repeats from the bases, where the catalog reads them from a file, and everything
/// after detection is deliberately the same code, so a disagreement points at the file
/// rather than at the policy.
///
/// # The whole contig is always classified, whatever `wanted` says (spec §2.5)
///
/// Classification is not local: a repeat at 999 is bundled away by a neighbour at 1030, so
/// a run that never looked past 1000 would call a locus the whole-contig run rejects —
/// same reference, different answers, because of `--regions`. So the five steps run over
/// the contig entire and `wanted` narrows only what comes out, which is
/// [`emit_into`]'s job.
///
/// # A rejected repeat is generic territory, not a hole (spec §2.2)
///
/// The generic path is the **default**; the other three are exceptions carved out
/// of it. So a repeat classification turns down for being impure, or low-copy, or
/// compound simply stays `Generic` — it is not a bundle, and it is certainly not a
/// hole. Only the *flank test* makes a bundle: a repeat with another repeat within
/// `bundle_threshold` of it, which is exactly the set [`segment_criteria::classify`] hands back.
pub fn partition_resident_in(
    chrom: &str,
    contig: ContigId,
    bases: &[u8],
    config: &TypedRegionConfig,
    wanted: &[GenomeRegion],
) -> (Vec<TypedRegion>, TypedRegionCounts) {
    let contig_len = bases.len() as u64;
    if contig_len == 0 {
        // A zero-length contig has no 1-based position to cover, so it has no
        // regions. `GenomeRegions` drops these before a caller gets here (C3); this
        // keeps the function total for a direct caller.
        return (Vec::new(), TypedRegionCounts::default());
    }
    let wanted: Vec<GenomeRegion> = wanted
        .iter()
        .filter(|span| span.contig == contig)
        .copied()
        .collect();

    // 1. Detect.
    let raw = find_tandem_repeats(bases, config.criteria.periods, &config.scan);
    // 2. Clean.
    let cleaned = segment_criteria::prefilter(&raw, &config.criteria);
    // 3. Admit — the whole contig, for the reason the doc gives.
    let classified = segment_criteria::classify(
        cleaned.clone(),
        chrom,
        bases,
        Position(1),
        Bp(contig_len),
        &config.criteria,
    );
    // 4. Cap: coverage runs over the *cleaned* intervals, then the satellite test.
    //
    //    Over the cleaned set, not the raw one (spec §2.4, §8): the raw scanner is
    //    deliberately permissive, so capping its coverage would let detector noise
    //    declare a satellite and silently swallow the real loci underneath it.
    //
    //    **Tested since D3** (`the_satellite_cap_applies_to_the_cleaned_coverage_not_
    //    the_raw`) — this comment used to say "known untested", and D1 was right that
    //    nothing then could discriminate the two. D3's 6 kb windowing fixture can: it
    //    carries real scanner noise (1181 raw runs vs 5 cleaned) and one noise interval
    //    abuts the array, so capping the raw set moves the satellite's edge.
    //
    //    D1's prediction of what a discriminating fixture needed — "≥ 1 kb of
    //    contiguous low-copy noise abutting a sub-cap array" — was too strong. Noise
    //    merely *touching* an over-cap run is enough to show the choice; the ≥ 1 kb case
    //    is what would show the *stakes* (noise inventing a satellite and swallowing the
    //    loci under it), and no fixture produces that yet.
    let runs = coverage_runs(&cleaned);

    // The tally's raw material, taken before the loci are consumed.
    //
    // **The repeat coverage is merged before it is measured**, because the count is of
    // bases and not of tracts: the cleaned set still holds intersecting tracts — the
    // pre-screen removes period-multiple re-detections of one tract, not two different
    // repeats that overlap — so summing the spans directly would charge a shared base once
    // per tract covering it. The satellite cap in this same function merges before it
    // measures ([`coverage_runs`]), and so does the catalog's tally.
    // Every field named, no `..default()` tail: a counter added later must break this line
    // rather than start life silently at zero in the catalog's own comparison.
    let mut counts = TypedRegionCounts {
        spans: wanted.len() as u64,
        ssr_loci: 0,
        ssr_bundles: 0,
        ssr_bundle_bp: 0,
        generic: 0,
        satellites: 0,
        satellite_bp: 0,
        repeat_bp_with_no_locus: covered_bp(&clipped_to_requested(
            &cleaned
                .iter()
                .map(|interval| (interval.start, interval.end))
                .collect::<Vec<_>>(),
            contig,
            &wanted,
        )),
        rejected_by_reason: RejectionCounts::default(),
    };
    // **A rejection is charged where the repeat starts**, and only if the caller asked
    // about that base — so a repeat is charged once however many spans can see it, and
    // the catalog charges the same way.
    for (interval, reason) in &classified.rejected {
        if inside_requested(interval.start, contig, &wanted) {
            counts
                .rejected_by_reason
                .add(*reason, interval.end - interval.start);
        }
    }

    // 5. Partition — the whole contig as **one block** (below), then fill every gap
    //    with `Generic`.
    let features = resolve_features(&runs, classified.loci, &classified.bundled, contig, config);

    let mut out = Vec::new();
    let mut emitted_upto = 0;
    for region in fill_generic_gaps(features, contig, contig_len) {
        emit_into(&mut out, region, &wanted, &mut emitted_upto);
    }
    for region in &out {
        charge_region(&mut counts, region, &wanted);
    }
    (out, counts)
}

/// Whether a region emitted past the edge of a requested span is **cut back to it**.
///
/// Only `Generic` is. A locus, a bundle and a satellite are each a claim about their own
/// extent, and half of one is a different claim: half a locus is not a locus, a clipped
/// bundle's members describe bases outside their own region, and a satellite clipped to
/// 100 bases contradicts the very cap that made it one. A generic stretch is the only kind
/// that is not a finding — *nothing more specific can be said here* stays true of any part
/// of it.
fn clips_at_a_bed_edge(kind: &RegionKind) -> bool {
    match kind {
        RegionKind::Generic => true,
        RegionKind::SsrSegment(_) | RegionKind::SsrBundle { .. } | RegionKind::Satellite => false,
    }
}

/// Narrow one region to what the caller asked for, and keep what survives (spec §2.5).
///
/// - **Outside every requested span** → dropped. It was classified so that the regions
///   inside would be *right*, not to be shown.
/// - **A finding** — locus, bundle or satellite ([`clips_at_a_bed_edge`]) → emitted
///   **whole**, even past the edge: the requested span grows to hold it, and that grown
///   span — the *effective* region — is what the partition invariant is stated over.
/// - **`Generic`** → **clipped** to each requested span it overlaps, which may be more
///   than one: every span on a contig is cut from that contig's one classification, so a
///   generic run across two of them must come back as two regions with the gap dropped,
///   not one region covering ground the caller did not ask for.
///
/// `emitted_upto` keeps the output non-overlapping when a finding has just been emitted
/// whole past an edge and the next requested span starts inside it.
fn emit_into(
    out: &mut Vec<TypedRegion>,
    region: TypedRegion,
    requested: &[GenomeRegion],
    emitted_upto: &mut u64,
) {
    let overlaps = |r: &GenomeRegion| {
        r.contig == region.region.contig
            && r.start.get() <= region.region.end.get()
            && region.region.start.get() <= r.end.get()
    };

    if !clips_at_a_bed_edge(&region.kind) {
        if requested.iter().any(overlaps) {
            *emitted_upto = region.region.end.get();
            out.push(region);
        }
        return;
    }

    for r in requested.iter().filter(|r| overlaps(r)) {
        let start = region
            .region
            .start
            .get()
            .max(r.start.get())
            .max(*emitted_upto + 1);
        let end = region.region.end.get().min(r.end.get());
        if start > end {
            continue;
        }
        *emitted_upto = end;
        out.push(TypedRegion {
            region: GenomeRegion {
                contig: region.region.contig,
                start: Position(start),
                end: Position(end),
            },
            kind: region.kind.clone(),
        });
    }
}

/// Charge one emitted region to the running tally.
fn charge_region(counts: &mut TypedRegionCounts, region: &TypedRegion, requested: &[GenomeRegion]) {
    let bp = region.region.len();
    match &region.kind {
        RegionKind::SsrSegment(_) => {
            counts.ssr_loci += 1;
            // This repeat coverage DID yield a locus, so it is not part of the gap.
            //
            // **Only the part inside the request is cancelled**, because only that part
            // was ever charged. A locus is emitted whole where it crosses a requested
            // edge, so subtracting its whole length would take back bases that were never
            // counted, and the counter would underflow at exactly the loci that straddle
            // one. What is inside the request is a subset of the coverage charged for it,
            // so this cannot.
            let inside: u64 = clipped_to_requested(
                &[(region.region.start.get() - 1, region.region.end.get())],
                region.region.contig,
                requested,
            )
            .into_iter()
            .map(|(start, end)| end - start)
            .sum();
            debug_assert!(inside <= bp, "a clipped locus cannot grow");
            debug_assert!(
                inside <= counts.repeat_bp_with_no_locus,
                "a locus cancels coverage that was charged for it, so this cannot underflow"
            );
            counts.repeat_bp_with_no_locus -= inside;
        }
        RegionKind::SsrBundle { .. } => {
            counts.ssr_bundles += 1;
            // **The number spec §10's bundle question needs and has never had** —
            // production drops these uncounted.
            counts.ssr_bundle_bp += bp;
        }
        RegionKind::Generic => counts.generic += 1,
        RegionKind::Satellite => {
            counts.satellites += 1;
            counts.satellite_bp += bp;
        }
    }
}

/// Step 5 for one **block**: cap the runs, place the loci the surviving runs do not
/// swallow, cluster the bundle members — the non-generic features, coordinate-ordered.
/// Generic is not this function's business: it is whatever is left over, and only the
/// caller knows how far "left over" reaches.
///
/// # How far a decision here can reach
///
/// Every rule this function applies has a **radius**: runs merge only when they abut
/// (radius 0), clustering chains members within `bundle_threshold` (that radius), and
/// swallowing is containment (radius 0). So a stretch of repeat structure separated from
/// the next by more than `bundle_threshold` of repeat-free sequence decides its own
/// contents and nothing else's, whichever way the caller arrived at the intervals.
///
/// That is what lets the catalog's reader hand this function a contig's stored rows and
/// get the answer a scan of the same contig gives.
///
/// **Both paths bottom out here, so the differential does not test it.** A bug in this
/// function is invisible to a comparison that runs it twice, by construction. It is
/// covered instead by the partition invariant and by `.cat` parity against the
/// trf-mod-built golden catalog (`repeat_catalog/anchor.rs`).
///
/// # Absorption: a satellite swallows what it touches (spec §2.4a)
///
/// **The rule (owner, 2026-07-16): a microsatellite or a cluster too close to a satellite
/// is swallowed by the satellite, which expands to cover it.**
///
/// A microsatellite 20 bp from a 1 kb array is not genotypeable, and the reason is the
/// array: there is no clean flank on that side, because the flank *is* array. A satellite
/// is already the region type that says *"an array — do not look for loci in here"* (spec
/// §2.1), so it is the one that should say it here too. The alternative — exempting the
/// array from the flank test so the neighbour becomes a clean locus — is wrong for a
/// right-sounding reason (spec §2.4).
///
/// **Both kinds, because both arise, and by different routes.** A cluster reaches here
/// when the array's own tract passes classification's gates and bundles with its neighbour. A
/// bare locus reaches here when it does *not*: a satellite run is built from **cleaned
/// coverage**, while bundling only ever sees tracts that cleared the scope/score/compound
/// gates — so an array rejected by one of those still forms a satellite and never bundles
/// anything. One rule covers both.
///
/// Absorption is iterated to a **fixed point** (a grown satellite reaches further than
/// the run it came from) and it subsumes spec §2.1's swallow: containment and adjacency
/// are the same predicate ([`absorb_into`]). Bundle members are absorbed **before** they
/// are clustered — an ordering that turned out to be load-bearing twice, once for the
/// rule and once for windowing; the reasons are at the loop.
///
/// **What this replaces, and why it was not a nicety.** The old test read a cluster's
/// **start** only, so the answer depended on which *side* of the array the microsat sat:
/// on the left it emitted a bundle **overlapping** the satellite — an invalid partition;
/// on the right the hull's start fell inside the run, so the cluster was dropped whole
/// and the microsat's bases silently became `Generic`. Same situation, two different
/// wrong answers. Probed before believing it, then fixed —
/// `a_microsatellite_beside_a_satellite_is_absorbed_into_it`.
pub(crate) fn resolve_features(
    runs: &[CoverageRun],
    loci: Vec<SsrSegment>,
    bundled: &[RepeatInterval],
    contig: ContigId,
    config: &TypedRegionConfig,
) -> Vec<TypedRegion> {
    let max_str_len = config.max_str_len.get();
    let bundle_threshold = config.criteria.bundle_threshold;
    let mut features: Vec<TypedRegion> = Vec::new();

    // The satellites: over-cap coverage runs — as spans that can still GROW (below).
    let mut satellites: Vec<CoverageRun> = runs
        .iter()
        .copied()
        .filter(|r| r.len() > max_str_len)
        .collect();

    let mut loci: Vec<(CoverageRun, SsrSegment)> = loci
        .into_iter()
        .map(|l| {
            (
                CoverageRun {
                    start: l.start(),
                    end: l.end(),
                },
                l,
            )
        })
        .collect();
    // 0-based half-open → 1-based inclusive, the same one conversion as everywhere
    // else (spec §4).
    let mut members: Vec<RepeatInterval> = bundled.to_vec();

    // **Absorption** (spec §2.4a; see the fn docs). Anything a satellite overlaps *or*
    // comes within `bundle_threshold` of is absorbed into it, and the satellite grows to cover
    // it. Iterated to a fixed point: a satellite that has grown reaches `bundle_threshold`
    // further than the run it came from, so it can absorb something the run could not.
    //
    // It terminates: each pass either absorbs at least one of a finite set of features
    // or stops. And it cannot reach out of the block — absorption's reach is the bundle
    // radius, which is narrower than the repeat-free sequence that bounds a block.
    //
    // **Bundle members are absorbed one by one, BEFORE they are clustered, and that
    // ordering is load-bearing twice over:**
    //
    // 1. It absorbs whole clusters, never part of one — so the survivors are still
    //    complete clusters and `bundle_clusters` can regroup them (its precondition:
    //    `bundled` is the concatenation of clusters in coordinate order). If a member
    //    is absorbed, every member chained to it is too: `is_close` implies "within
    //    `bundle_threshold`" (each of its four clauses puts a pair of endpoints inside the
    //    radius), so the grown satellite reaches the next member in the chain, and so
    //    on by induction.
    // 2. It is what makes the **windowed** walk correct at all. A window truncates a
    //    detection at its scanned slice's edge, and only a tract longer than
    //    `max_str_len` can be truncated — i.e. an array. A truncated member's `end`
    //    is wrong, so `is_close` cannot re-chain it, and `bundle_clusters` sees a
    //    singleton. Absorbing members first means a truncated one never reaches
    //    clustering: it is over-cap, so it lies inside its own satellite run.
    //    (`bundle_clusters`' `debug_assert` is what caught this, exactly as D2 built it
    //    to.)
    loop {
        let mut absorbed = false;
        members.retain(|iv| {
            let span = CoverageRun {
                start: iv.start + 1,
                end: iv.end,
            };
            !absorb_into(&mut satellites, span, bundle_threshold, &mut absorbed)
        });
        loci.retain(|(span, _)| {
            !absorb_into(&mut satellites, *span, bundle_threshold, &mut absorbed)
        });
        if !absorbed {
            break;
        }
    }

    // Bundles (D2). `classify` set these aside — repeats too close to each other for
    // any of them to have a clean flank (spec §2.4) — and handed them back rather
    // than deleting them, which is the whole point of `Classified::bundled`. Each
    // cluster becomes **one** region spanning the hull of its tracts: the gaps
    // between members are inside it, and rightly so — they are shorter than a
    // flank, so nothing can be anchored in them either.
    //
    // A surviving cluster cannot touch a satellite even though only its *members* were
    // tested: its hull adds only the gaps between members, and those are shorter than a
    // flank — far too short to hide an over-cap run — so a hull that reaches a satellite
    // has a member that reaches it.
    let clusters: Vec<(CoverageRun, Vec<RepeatInterval>)> =
        segment_criteria::bundle_clusters(&members, bundle_threshold)
            .into_iter()
            .map(|cluster| {
                let hull = CoverageRun {
                    start: cluster.first().expect("non-empty cluster").start + 1,
                    end: cluster.iter().map(|iv| iv.end).max().expect("non-empty"),
                };
                (hull, cluster)
            })
            .collect();

    for satellite in &satellites {
        features.push(TypedRegion {
            region: GenomeRegion {
                contig,
                start: Position(satellite.start),
                end: Position(satellite.end),
            },
            kind: RegionKind::Satellite,
        });
    }
    for (span, locus) in loci {
        features.push(TypedRegion {
            region: GenomeRegion {
                contig,
                start: Position(span.start),
                end: Position(span.end),
            },
            kind: RegionKind::SsrSegment(locus),
        });
    }
    for (hull, cluster) in clusters {
        features.push(TypedRegion {
            region: GenomeRegion {
                contig,
                start: Position(hull.start),
                end: Position(hull.end),
            },
            kind: RegionKind::SsrBundle {
                tracts: cluster.into_boxed_slice(),
            },
        });
    }

    features.sort_by_key(|f| f.region.start);
    features
}

/// Absorb `span` into any satellite it overlaps or lies within `bundle_threshold` of: those
/// satellites and `span` become **one** span covering all of them, gaps included.
/// Returns whether it was absorbed (and sets `absorbed`, the fixed-point loop's flag).
///
/// # Why the gap, and not `segment_criteria::is_close`
///
/// This asks a **flank** question — *are there `bundle_threshold` clean bases between this
/// feature and the array?* — so the gap is the measure, and `< bundle_threshold` is strict, as
/// classification's own gap clause is (`is_close_is_strict_at_the_threshold`).
///
/// `is_close` itself is the wrong tool here despite testing the same relation between
/// two *tracts*: it is four `abs_diff` clauses ported from GangSTR, and three of them
/// compare start-to-start and end-to-end — meaningful between two 25 bp tracts, noise
/// against a span that can be 2 Mb long, where a locus 20 bp past the end is millions
/// of bases from the start. The relation is the same; the predicate cannot be.
///
/// The gaps swept in are deliberate: fewer than `bundle_threshold` bases between two repeats is
/// sequence nothing can be anchored in — the same reasoning that puts a cluster's
/// internal gaps inside its hull (spec §2.4).
fn absorb_into(
    satellites: &mut Vec<CoverageRun>,
    span: CoverageRun,
    bundle_threshold: u64,
    absorbed: &mut bool,
) -> bool {
    // `span.start <= s.end + bundle_threshold` is "the gap on this side is < bundle_threshold", and it
    // is also true whenever the two overlap — so containment (a locus *inside* a
    // satellite, spec §2.1's swallow) and adjacency are one rule, not two. They were
    // two, and reading only the hull's start is what made the answer depend on which
    // side of the array the feature sat.
    let touching: Vec<usize> = satellites
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            span.start <= s.end + bundle_threshold && s.start <= span.end + bundle_threshold
        })
        .map(|(i, _)| i)
        .collect();
    if touching.is_empty() {
        return false;
    }

    let mut union = span;
    // Back to front: removing by index cannot disturb an index still to come.
    for &i in touching.iter().rev() {
        let s = satellites.remove(i);
        union.start = union.start.min(s.start);
        union.end = union.end.max(s.end);
    }
    satellites.push(union);
    merge_runs(satellites);
    *absorbed = true;
    true
}

/// Union overlapping and abutting runs in place, leaving them ascending.
///
/// Shared by the satellite absorption above and by [`coverage_runs`] — the same rule
/// (touching spans are one span), so the same code.
fn merge_runs(runs: &mut Vec<CoverageRun>) {
    runs.sort_by_key(|s| (s.start, s.end));
    let mut merged: Vec<CoverageRun> = Vec::with_capacity(runs.len());
    for s in runs.drain(..) {
        match merged.last_mut() {
            // `s.start <= last.end + 1` merges abutting runs as well as overlapping ones.
            Some(last) if s.start <= last.end + 1 => last.end = last.end.max(s.end),
            _ => merged.push(s),
        }
    }
    *runs = merged;
}

/// How many **distinct** bases a set of half-open spans covers, a base shared by two spans
/// counting once.
///
/// The spans arrive in the detector's order — period-major, not coordinate order — and may
/// intersect, so they are sorted and unioned rather than summed.
fn covered_bp(spans: &[(u64, u64)]) -> u64 {
    let mut spans = spans.to_vec();
    spans.sort_unstable();
    let mut total = 0;
    let mut open: Option<(u64, u64)> = None;
    for (start, end) in spans {
        match open {
            // Sorted by start, so `start >= run_start` and the union is just the wider end.
            Some((run_start, run_end)) if start <= run_end => {
                open = Some((run_start, run_end.max(end)));
            }
            Some((run_start, run_end)) => {
                total += run_end - run_start;
                open = Some((start, end));
            }
            None => open = Some((start, end)),
        }
    }
    if let Some((run_start, run_end)) = open {
        total += run_end - run_start;
    }
    total
}

/// The parts of `spans` (0-based half-open, on `contig`) that fall inside `wanted`.
///
/// `wanted` is 1-based inclusive, ng's convention for a requested region, so each is
/// `[start - 1, end)` in the detector's space.
fn clipped_to_requested(
    spans: &[(u64, u64)],
    contig: ContigId,
    wanted: &[GenomeRegion],
) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    for &(start, end) in spans {
        for region in wanted.iter().filter(|r| r.contig == contig) {
            let from = start.max(region.start.get() - 1);
            let to = end.min(region.end.get());
            if from < to {
                out.push((from, to));
            }
        }
    }
    out
}

/// Whether `position` (0-based) lies inside one of `wanted` on `contig`.
fn inside_requested(position: u64, contig: ContigId, wanted: &[GenomeRegion]) -> bool {
    wanted
        .iter()
        .any(|r| r.contig == contig && position >= r.start.get() - 1 && position < r.end.get())
}

/// A merged run of repeat coverage, 1-based inclusive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct CoverageRun {
    start: u64,
    end: u64,
}

impl CoverageRun {
    /// The run's first base, 1-based inclusive.
    ///
    /// Shared with the catalog for the same reason [`Self::len`] is: a region-scoped tally
    /// clips these runs to the spans the caller asked for, and both sides must clip the same
    /// runs the same way.
    pub(crate) fn start(self) -> u64 {
        self.start
    }

    /// The run's last base, 1-based inclusive.
    pub(crate) fn end(self) -> u64 {
        self.end
    }

    /// The run's length in bases, inclusive at both ends.
    ///
    /// **Shared with the catalog's derived segmentation**
    /// ([`crate::ng::repeat_catalog::segments`]), which sums it over a contig's runs to get
    /// the repeat coverage charged to `repeat_bp_with_no_locus`. Same rule, same
    /// arithmetic, so the two tallies cannot drift on it.
    pub(crate) fn len(self) -> u64 {
        self.end - self.start + 1
    }
}

/// Union the intervals into maximal runs of covered bases, 1-based inclusive.
///
/// Overlapping **and abutting** runs merge: two tracts that touch cover a
/// contiguous stretch of repeat, and whether that stretch is a satellite is a
/// question about the stretch, not about how the detector happened to split it.
///
/// Input is `RepeatInterval`'s 0-based half-open; output is ng's 1-based
/// inclusive, so `[s, e)` becomes `[s + 1, e]` — the same one conversion `classify`
/// makes (spec §4).
pub(crate) fn coverage_runs(intervals: &[RepeatInterval]) -> Vec<CoverageRun> {
    let mut spans: Vec<CoverageRun> = intervals
        .iter()
        .filter(|iv| iv.end > iv.start)
        .map(|iv| CoverageRun {
            start: iv.start + 1,
            end: iv.end,
        })
        .collect();
    merge_runs(&mut spans);
    spans
}

/// Fill every gap between `features` with `Generic`, so the result tiles
/// `[1, contig_len]` exactly.
///
/// **Maximality is a correctness requirement here, not tidiness** (spec §2.3): a
/// generic region is territory the pileup mints loci *inside*, so its reach is
/// bounded by the region it was handed. Split a run at *p* and an indel spanning
/// *p* is callable by neither half — it never appears, and nothing fails. Hence
/// one `Generic` per gap, however long.
///
/// `features` must be coordinate-ordered and non-overlapping.
pub(crate) fn fill_generic_gaps(
    features: Vec<TypedRegion>,
    contig: ContigId,
    contig_len: u64,
) -> Vec<TypedRegion> {
    let generic = |start: u64, end: u64| TypedRegion {
        region: GenomeRegion {
            contig,
            start: Position(start),
            end: Position(end),
        },
        kind: RegionKind::Generic,
    };

    let mut out = Vec::with_capacity(features.len() * 2 + 1);
    let mut pos = 1u64;
    for f in features {
        if f.region.start.get() > pos {
            out.push(generic(pos, f.region.start.get() - 1));
        }
        pos = f.region.end.get() + 1;
        out.push(f);
    }
    if pos <= contig_len {
        out.push(generic(pos, contig_len));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use segment_criteria::DEFAULT_BUNDLE_THRESHOLD;

    /// **A base under two repeats is one base.** The pre-screen leaves intersecting tracts
    /// in place — it removes period-multiple re-detections of *one* tract, not two different
    /// repeats that happen to overlap — so the accumulator behind
    /// `repeat_bp_with_no_locus` has to union before it measures. It summed instead until
    /// the catalog's tally was compared against it and came out a base short over 200 kb.
    ///
    /// The spans arrive period-major rather than in coordinate order, which is why the last
    /// case is out of order on purpose.
    #[test]
    fn overlapping_coverage_counts_a_shared_base_once() {
        assert_eq!(covered_bp(&[]), 0);
        assert_eq!(covered_bp(&[(10, 20)]), 10);
        // Disjoint: 10 + 5.
        assert_eq!(covered_bp(&[(10, 20), (30, 35)]), 15);
        // Overlapping by one base: 10 + 10 - 1, not 20.
        assert_eq!(covered_bp(&[(10, 20), (19, 29)]), 19);
        // Abutting is not overlapping: 10 + 10.
        assert_eq!(covered_bp(&[(10, 20), (20, 30)]), 20);
        // One span swallowed by another, and handed over in the wrong order.
        assert_eq!(covered_bp(&[(12, 15), (10, 30)]), 20);
    }

    /// `Default` is the **short-read** settings now (spec §2.3), not the catalog's —
    /// which is what the golden oracles pin explicitly (see [`catalog_config`]).
    /// Pinned against the **literals**, not against the consts that define it:
    /// asserting `periods.min() == DEFAULT_MIN_PERIOD` would compare a constant with
    /// itself and could not fail (the same tautology the A2 review caught).
    #[test]
    fn default_config_is_the_short_read_settings() {
        let c = TypedRegionConfig::default();
        // One period range (scanner detects = classification accepts), living in
        // `criteria`. §2.3 lowered the floor to 1: a qualifying homopolymer is a
        // period-1 STR locus.
        assert_eq!(
            c.criteria.periods.min(),
            1,
            "period-1 homopolymers classified"
        );
        assert_eq!(c.criteria.periods.max(), 6, "..to the hexa ceiling");
        assert_eq!(c.max_str_len, Bp(100), "the short-read satellite cap");
        assert_eq!(c.scan, ScanParams::default());
        assert_eq!(c.criteria, SsrSegmentCriteria::default());
        // The short-read copy-number floors `[8,6,6,6,5,4]` (spec ssr §5.1): the mono floor
        // drops to 6 (Illumina read-artifact onset), the di floor to 4.
        let floors: Vec<u32> = (1..=6)
            .map(|p| c.criteria.min_copies.for_period(p))
            .collect();
        assert_eq!(floors, vec![8, 6, 6, 6, 5, 4], "the short-read floors");
    }

    /// The catalog's settings — what [`TypedRegionConfig::default`] carried before
    /// spec §2.3 moved `Default` to the short-read floors. The `.cat` parity oracle
    /// and the fixtures below were written against these (di..hexa,
    /// copy floors `[10,5,4,3,3,3]`, a **50 bp** flank/bundle radius, a **1 kb**
    /// satellite cap), so they pin them explicitly now that `Default` no longer
    /// supplies them (spec §8.1: the oracle must be pinned to the catalog's
    /// settings, not to whatever `Default` is).
    fn catalog_config() -> TypedRegionConfig {
        use crate::ng::tandem_repeat::PeriodRange;
        use segment_criteria::MinCopies;
        TypedRegionConfig {
            max_str_len: Bp(1000),
            criteria: SsrSegmentCriteria {
                periods: PeriodRange::new(2, 6).expect("2..=6 is a valid period range"),
                min_copies: MinCopies::new([10, 5, 4, 3, 3, 3], 3),
                bundle_threshold: 50,
                ..SsrSegmentCriteria::default()
            },
            ..TypedRegionConfig::default()
        }
    }

    #[test]
    fn counts_start_at_zero() {
        let c = TypedRegionCounts::default();
        assert_eq!(c.spans, 0);
        assert_eq!(c.ssr_loci, 0);
        assert_eq!(c.ssr_bundle_bp, 0);
        assert_eq!(c.repeat_bp_with_no_locus, 0);
    }

    // ---- GenomeRegions (C3) ---------------------------------------------

    const CONTIGS: &[ContigBounds] = &[
        ContigBounds {
            name: "chr1",
            length: 100,
        },
        ContigBounds {
            name: "chr2",
            length: 50,
        },
    ];

    /// `whole_contigs` is the default, and the spans are what `regions.rs` calls
    /// "the region set whose every region covers an entire contig" — full-length,
    /// **1-based inclusive**, one per contig, in table order.
    #[test]
    fn whole_contigs_covers_each_contig_end_to_end() {
        let g = GenomeRegions::whole_contigs(CONTIGS);
        let regions: Vec<_> = g.iter().collect();

        assert_eq!(g.len(), 2);
        assert!(!g.is_empty());
        assert_eq!(regions[0].contig, ContigId(0));
        assert_eq!(regions[0].start, Position(1), "1-based: starts at 1, not 0");
        assert_eq!(
            regions[0].end,
            Position(100),
            "inclusive: the last base IS 100"
        );
        assert_eq!(regions[0].len(), 100, "a 100 bp contig walks 100 bases");
        assert_eq!(regions[1].contig, ContigId(1));
        assert_eq!(regions[1].end, Position(50));
        assert_eq!(regions[1].len(), 50);
    }

    /// **No rebasing happens here, and that is the finding.** Spec §4 expected this
    /// seam to widen *and rebase*; `regions::Region` is already 1-based inclusive
    /// (its own invariant is `1 <= start <= end`), so only the width converts.
    ///
    /// This test is the guard on that: if production's base ever moved, the
    /// coordinates below would shift by one and ng's whole 1-based contract would
    /// quietly break at its busiest boundary.
    #[test]
    fn the_seam_widens_but_does_not_rebase() {
        let production = RegionSet::whole_contigs(CONTIGS);
        let ours = GenomeRegions::whole_contigs(CONTIGS);

        for (p, n) in production.iter().zip(ours.iter()) {
            assert_eq!(
                u64::from(p.start),
                n.start.get(),
                "start is carried across verbatim — production is already 1-based"
            );
            assert_eq!(u64::from(p.end), n.end.get(), "end likewise");
            assert_eq!(p.chrom_id, n.contig.get(), "and the id is the same index");
        }
    }

    /// Zero-length contigs contribute no span, so they never reach a consumer —
    /// `RegionSet`'s rule, inherited. This is why spec §2.3 can assert "zero-length
    /// contigs are never asked about" without a consumer guarding for it.
    #[test]
    fn a_zero_length_contig_is_dropped_before_a_consumer_sees_it() {
        let contigs = &[
            ContigBounds {
                name: "empty",
                length: 0,
            },
            ContigBounds {
                name: "chr1",
                length: 10,
            },
        ];
        let g = GenomeRegions::whole_contigs(contigs);
        let regions: Vec<_> = g.iter().collect();
        assert_eq!(regions.len(), 1, "the empty contig contributes nothing");
        assert_eq!(
            regions[0].contig,
            ContigId(1),
            "and the ids do NOT renumber"
        );
    }

    /// A BED round-trip: ng inherits `RegionSet`'s parsing, its 0-based-BED → 1-based
    /// conversion, and its coalescing — none of which ng reimplements. The
    /// overlapping pair must come back merged.
    #[test]
    fn from_bed_path_parses_converts_and_coalesces() {
        use std::io::Write;
        std::fs::create_dir_all("tmp").unwrap();
        let dir = tempfile::tempdir_in("tmp").unwrap();
        let bed = dir.path().join("r.bed");
        {
            let mut f = std::fs::File::create(&bed).unwrap();
            // BED is 0-based half-open: [0,10) is 1-based [1,10].
            writeln!(f, "chr1\t0\t10").unwrap();
            // Overlaps the first — must coalesce into [1, 20].
            writeln!(f, "chr1\t5\t20").unwrap();
            writeln!(f, "chr2\t0\t5").unwrap();
        }
        let g = GenomeRegions::from_bed_path(&bed, CONTIGS).expect("valid bed");
        let regions: Vec<_> = g.iter().collect();

        assert_eq!(regions.len(), 2, "the chr1 pair coalesced");
        assert_eq!(regions[0].contig, ContigId(0));
        assert_eq!(
            (regions[0].start, regions[0].end),
            (Position(1), Position(20)),
            "BED 0-based [0,20) becomes 1-based inclusive [1,20]"
        );
        assert_eq!(regions[1].contig, ContigId(1));
        assert_eq!(
            (regions[1].start, regions[1].end),
            (Position(1), Position(5))
        );
    }

    // ---- The resident partition ------------------------------------------

    /// **The invariant — the acceptance test** (spec §2.3).
    ///
    /// Within a requested region the typed regions are **contiguous**
    /// (`start == prev.end + 1`), **non-overlapping**, **complete** (their union is
    /// the whole span), and **maximal** (no two consecutive share a kind).
    ///
    /// One property: *concatenating the regions reconstructs what was asked for,
    /// exactly.* Every way this design fails shows up as a violation — a rejected
    /// repeat left as a hole breaks completeness; a flank counted as ownership
    /// breaks non-overlap; a generic run split at a window edge breaks maximality.
    #[track_caller]
    fn assert_partitions(regions: &[TypedRegion], contig: ContigId, contig_len: u64, case: &str) {
        assert!(
            !regions.is_empty(),
            "{case}: a non-empty contig has regions"
        );
        let mut expected_start = 1u64;
        let mut prev_kind: Option<std::mem::Discriminant<RegionKind>> = None;
        for r in regions {
            assert_eq!(r.region.contig, contig, "{case}: contig");
            assert_eq!(
                r.region.start.get(),
                expected_start,
                "{case}: gap or overlap at {} (expected {expected_start}); regions: {regions:#?}",
                r.region.start.get()
            );
            assert!(
                r.region.end >= r.region.start,
                "{case}: empty region {:?}",
                r.region
            );
            let kind = std::mem::discriminant(&r.kind);
            assert_ne!(
                Some(kind),
                prev_kind,
                "{case}: two consecutive regions share a kind at {} — MAXIMALITY. \
                 For Generic this is a correctness bug, not untidiness: the pileup \
                 mints loci inside a Generic region, so a split run makes an indel \
                 across the join callable by neither half.",
                r.region.start.get()
            );
            prev_kind = Some(kind);
            expected_start = r.region.end.get() + 1;
        }
        assert_eq!(
            expected_start - 1,
            contig_len,
            "{case}: the partition must cover exactly [1, {contig_len}] — COMPLETENESS"
        );
    }

    /// A contig with one clean isolated (AT)*8 tract: Generic / SsrSegment / Generic.
    /// The smallest case that shows the partition doing its job.
    ///
    /// **The flanks must be aperiodic, and the first version's were not** (fixed
    /// 2026-07-20). They were `(CGCA)*15` — a period-4 tandem repeat — so the contig was
    /// not "one lone tract" at all: the scanner read the whole thing as a single impure
    /// period-4 tract (`0..136 p4`, the 16 bp `AT` insert being exactly four periods of
    /// four) with the `AT` tract nested inside it. The assertion passed only because the
    /// then-current `is_close` could not see containment, so the nested `AT` came out a
    /// standalone locus and the enclosing p4 tract was dropped for want of a flank at the
    /// contig edge — two bugs cancelling. With [`segment_criteria::joins_cluster`] the
    /// nesting is seen and the honest answer for *that* fixture is one `SsrBundle`, which
    /// is how the fixture was caught. The flanks below are aperiodic, so the contig now
    /// holds the one tract this test is named for.
    #[test]
    fn a_lone_tract_partitions_as_generic_locus_generic() {
        // 60 bp of aperiodic sequence either side, so the default 50 bp flanks fit and
        // nothing in them is itself a tract (see the note above).
        let mut bases = b"ACGTTGCAAGCTCCTAGGATCGATTGCACGGTACCTGAAGCTTGCACTGATCCGTAGGCA".to_vec();
        let tract_start_0 = bases.len();
        bases.extend_from_slice(b"ATATATATATATATAT"); // 8 copies
        bases.extend_from_slice(b"TGCATTGGACCTAAGCGTTCAGGCTTACGATCCAGGTTACGATCCAAGTGCTTAGCATCG");
        let len = bases.len() as u64;

        let regions = partition_resident("chr1", ContigId(0), &bases, &catalog_config());
        assert_partitions(&regions, ContigId(0), len, "lone tract");

        let kinds: Vec<_> = regions
            .iter()
            .map(|r| match &r.kind {
                RegionKind::SsrSegment(_) => "locus",
                RegionKind::SsrBundle { .. } => "bundle",
                RegionKind::Generic => "generic",
                RegionKind::Satellite => "satellite",
            })
            .collect();
        assert_eq!(kinds, vec!["generic", "locus", "generic"]);

        // And the locus is where the tract is — 1-based, so 0-based + 1.
        let RegionKind::SsrSegment(l) = &regions[1].kind else {
            unreachable!()
        };
        assert_eq!(l.start(), tract_start_0 as u64 + 1);
        assert_eq!(l.motif().as_bytes(), b"AT");
        assert_eq!(regions[1].region.start.get(), l.start(), "region == tract");
    }

    /// A 2 kb array is **one** `Satellite`, and the locus inside it is **swallowed**
    /// — typed as one object, not searched for loci inside (spec §2.1).
    ///
    /// **The swallow has to be checked positively**, and this test's first version
    /// did not: it used `vec![b'C'; 60]` "flanks", which are a period-1 homopolymer
    /// — so the scanner found one period-2 tract spanning the *whole contig*, which
    /// starts at base 1, has no left flank, and was dropped. Nothing was swallowed
    /// because nothing was ever classified, and the assertion passed for entirely the
    /// wrong reason (mutation caught it: "don't drop loci inside a satellite"
    /// survived). The control below is the fix: at the same settings but a cap
    /// above the array, the same bases DO yield a locus.
    #[test]
    fn a_long_array_is_one_satellite_and_swallows_the_locus_inside_it() {
        let mut bases = filler(60);
        for _ in 0..1000 {
            bases.extend_from_slice(b"AT"); // 2 kb, well over the satellite cap
        }
        bases.extend(filler(60));
        let len = bases.len() as u64;
        let config = TypedRegionConfig::default();

        let regions = partition_resident("chr1", ContigId(0), &bases, &config);
        assert_partitions(&regions, ContigId(0), len, "satellite");

        let satellites: Vec<_> = regions
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::Satellite))
            .collect();
        assert_eq!(satellites.len(), 1, "one array, ONE satellite region");
        assert!(satellites[0].region.len() >= 2000);
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
            "the locus is swallowed by the satellite — the §2.4 cost, made visible"
        );

        // **The control.** Raise the cap above the array and the SAME bases classify a
        // locus — so the absence above is the cap doing its job, not classification
        // quietly rejecting the tract for some unrelated reason. This is also what
        // makes `max_str_len` a parameter rather than a fact of nature (spec §10).
        let uncapped = TypedRegionConfig {
            max_str_len: Bp(10_000),
            ..config
        };
        let regions = partition_resident("chr1", ContigId(0), &bases, &uncapped);
        assert_partitions(&regions, ContigId(0), len, "satellite, uncapped");
        assert_eq!(
            regions
                .iter()
                .filter(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
                .count(),
            1,
            "above the cap the very same tract IS a locus"
        );
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::Satellite)),
            "and nothing exceeds the raised cap, so no satellite"
        );
    }

    /// Coverage runs merge **abutting** tracts, not only overlapping ones: two
    /// tracts that touch cover a contiguous stretch of repeat, and whether *that*
    /// is a satellite is a question about the stretch, not about where the detector
    /// happened to split it.
    ///
    /// Untested until mutation said so — no other fixture puts two runs exactly
    /// end to end.
    #[test]
    fn abutting_coverage_runs_merge_into_one() {
        let iv = |start, end, period| RepeatInterval {
            start,
            end,
            period,
            score: 1,
        };
        // 1-based: [1,10] and [11,20] touch → one run [1,20]. [30,40] is separate.
        let runs = coverage_runs(&[iv(0, 10, 2), iv(10, 20, 3), iv(29, 40, 2)]);
        assert_eq!(
            runs,
            vec![
                CoverageRun { start: 1, end: 20 },
                CoverageRun { start: 30, end: 40 },
            ],
            "abutting runs merge; a gap of even one base does not"
        );
        assert_eq!(runs[0].len(), 20, "inclusive length");

        // The consequence, and the reason the merge exists: two 600 bp tracts that
        // touch are a 1200 bp run — a satellite — though neither tract is.
        let runs = coverage_runs(&[iv(0, 600, 2), iv(600, 1200, 2)]);
        assert_eq!(runs.len(), 1);
        assert_eq!(
            runs[0].len(),
            1200,
            "over the satellite cap, though neither tract alone is"
        );
    }

    /// A contig with a **3 kb array** — three times the 1 kb margin, so a walk that scanned
    /// only a grown BED span could not see the whole of it.
    fn contig_with_a_long_array() -> Vec<u8> {
        let mut bases = filler(8000);
        for i in 0..1500 {
            bases[1000 + i * 2..1002 + i * 2].copy_from_slice(b"AT");
        }
        bases
    }

    /// The kinds, in order — for asserting a partition's shape without spelling out four
    /// `matches!` arms at every site.
    fn kinds_of(regions: &[TypedRegion]) -> Vec<&'static str> {
        regions
            .iter()
            .map(|r| match &r.kind {
                RegionKind::SsrSegment(_) => "locus",
                RegionKind::SsrBundle { .. } => "bundle",
                RegionKind::Generic => "generic",
                RegionKind::Satellite => "satellite",
            })
            .collect()
    }

    /// Aperiodic filler — **not** a homopolymer, and that matters.
    ///
    /// A `vec![b'C'; 60]` "flank" is a period-1 tract, and worse: with an `(AT)n`
    /// array between two of them the scanner finds a **single period-2 tract
    /// spanning the whole contig**, which then starts at base 1, has no left flank,
    /// and is dropped — so a test asserting "no locus here" passes for entirely the
    /// wrong reason. This filler has no repeat at any period 1..=6 (see
    /// `a_repeat_free_contig_is_one_generic_region`, which is one `Generic` over
    /// it).
    fn filler(n: usize) -> Vec<u8> {
        b"ACGTTGCAAGCTTGCA"
            .iter()
            .copied()
            .cycle()
            .take(n)
            .collect()
    }

    /// A repeat-free contig is exactly one `Generic` region. Maximality: not many.
    #[test]
    fn a_repeat_free_contig_is_one_generic_region() {
        let bases = b"ACGTTGCAAGCTTGCAACGTTGCAAGCTTGCAACGTTGCAAGCTTGCA".repeat(3);
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());

        assert_partitions(&regions, ContigId(0), len, "repeat-free");
        assert_eq!(regions.len(), 1, "one Generic, not a run of them");
        assert!(matches!(regions[0].kind, RegionKind::Generic));
        assert_eq!(regions[0].region.start, Position(1));
        assert_eq!(regions[0].region.end, Position(len));
    }

    /// **An impure tract is `Generic`** (spec §8's fixture list, §2.2) — with the control
    /// that shows it is the interruption's doing.
    ///
    /// # But "impure → Generic" is not one rule, and it took this test to see it
    ///
    /// The scanner decides an impure tract's fate long before classification does, because
    /// Ruzzo–Tompa returns **maximal-scoring** segments. So an interruption has three
    /// possible outcomes, and only the third is this fixture:
    ///
    /// 1. **small** — the surrounding matches pay for it, the tract stays whole, and it is
    ///    pure enough to classify: a **locus**;
    /// 2. **large, with long pieces** — the segment splits, and the pure pieces are close
    ///    together: a **bundle** (or two loci, if far apart);
    /// 3. **large, with short pieces** — the segment splits and each piece falls under the
    ///    copy floor: **`Generic`**, which is this fixture (two 8 bp halves = 4 copies,
    ///    under the floor of 5).
    ///
    /// **Admission's purity gate is not what does any of this.** It is unreachable from the
    /// walk: a tract impure enough to fail the 0.80 floor always contains a purer
    /// sub-segment that scores higher, so Ruzzo–Tompa emits *that* instead — measured, not
    /// argued (a 0.79-purity fixture comes back a **locus**, its pure core). See
    /// `segment_criteria::the_walk_reaches_only_one_of_classifications_five_gates`.
    #[test]
    fn an_impure_tract_is_generic_when_its_pieces_fall_under_the_copy_floor() {
        let mut bases = filler(240);
        bases[100..126].copy_from_slice(b"ATATATATGGGGGGGGGGATATATAT");
        let len = bases.len() as u64;

        let regions = partition_resident("chr1", ContigId(0), &bases, &catalog_config());
        assert_partitions(&regions, ContigId(0), len, "impure tract");
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
            "an impure tract is not a locus: {:?}",
            kinds_of(&regions)
        );
        // Not a hole either — spec §2.2's property, which is what the fixture is for.
        assert_eq!(regions.len(), 1, "the whole contig is one Generic region");

        // **The control.** Remove the interruption and the same span at the same settings
        // IS a locus — so `Generic` above is the interruption's doing, not a tract that was
        // never admissible (D1's trap).
        let mut bases = filler(240);
        bases[100..126].copy_from_slice(b"ATATATATATATATATATATATATAT");
        let regions = partition_resident("chr1", ContigId(0), &bases, &catalog_config());
        assert_eq!(
            kinds_of(&regions),
            vec!["generic", "locus", "generic"],
            "pure, the same span is a locus"
        );
    }

    /// **A homopolymer is `Generic` at periods 2..=6** (spec §8's fixture list): nothing
    /// classifies period 1, and the bases are still covered. Pinned at the catalog's
    /// period scope explicitly via [`catalog_config`] — as of spec §2.3 ng's `Default`
    /// is `--min-period 1`, where a homopolymer IS a locus
    /// (`a_homopolymer_of_six_or_more_is_a_period_one_locus_at_default` pins that).
    ///
    /// `prefilter` is what removes it, and **the whole homopolymer, not just its period-1
    /// label** — the 2026-07-17 ordering fix. A poly-A tiles under `AA`, `AAA`, `AAAAA`, so
    /// the scanner emits the same span at every period in scope; only the period-1 interval
    /// divides them all, so it is kept as an eliminator through redundancy elimination and
    /// dropped by the period floor afterwards. Before that fix, periods 2, 3 and 5 survived
    /// and one homopolymer entered the partition as three "repeats"
    /// (`a_homopolymer_does_not_survive_as_a_period_two_repeat` pins it at the unit level).
    ///
    /// The neighbouring tract is the other half: the copy floor drops period-1 specks
    /// before they can eliminate anything, so a real STR is never taken out by one — the
    /// poly-A cascade.
    #[test]
    fn a_homopolymer_is_generic_and_does_not_take_its_neighbour_with_it() {
        let mut bases = filler(300);
        // 20 bp of poly-A, then a clean (CAG) tract 60 bp away — far enough not to bundle.
        bases[100..120].copy_from_slice(&[b'A'; 20]);
        bases[180..204].copy_from_slice(b"CAGCAGCAGCAGCAGCAGCAGCAG");
        let len = bases.len() as u64;

        let regions = partition_resident("chr1", ContigId(0), &bases, &catalog_config());
        assert_partitions(&regions, ContigId(0), len, "homopolymer");

        let loci: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) => Some(l),
                _ => None,
            })
            .collect();
        assert_eq!(
            loci.len(),
            1,
            "exactly one locus: the (CAG) tract. The homopolymer is not one, and — the \
             point — it did not eliminate the tract either: {:#?}",
            kinds_of(&regions)
        );
        assert_eq!(loci[0].motif().as_bytes(), b"CAG");
        assert!(
            loci[0].start() >= 181 && loci[0].start() <= 185,
            "and it is the tract at ~181, not something the cascade left behind: {}",
            loci[0].start()
        );
        // The homopolymer's own bases are covered, not a hole (spec §2.2).
        assert!(
            regions
                .iter()
                .any(|r| r.region.contains(Position(110)) && matches!(r.kind, RegionKind::Generic)),
            "the poly-A run is Generic territory"
        );

        // **And it reaches classification not at all** — it is gone by the pre-filter, so it is
        // turned down by no gate and counted under no reason. That is what "out of scope"
        // has to mean: before the 2026-07-17 fix this fixture recorded THREE `Compound`
        // rejections for one homopolymer, because its period-2/3/5 aliases each reached
        // `classify` separately. A rejection count that moves with the number of divisors of
        // a run's length is not measuring anything.
        let counts = counts_over_with(&bases, catalog_config());
        assert_eq!(
            counts.rejected_by_reason,
            RejectionCounts::default(),
            "the homopolymer is filtered out, not rejected — and the (CAG) tract is \
             classified, so nothing here is turned down at all: {:?}",
            counts.rejected_by_reason
        );
    }

    /// **At `Default` (spec §2.3) a homopolymer of ≥ 6 bp IS a period-1 STR locus** —
    /// the mirror of `a_homopolymer_is_generic_and_does_not_take_its_neighbour_with_it`.
    /// `--min-period 1` with a mononucleotide copy floor of 6 is the short-read
    /// default, so a poly-A run routes to STR handling instead of falling through as
    /// `Generic`. §2.3 flags this as the change that most alters the partition —
    /// genomes are dense in homopolymers, so expect many more `SsrSegment`s.
    #[test]
    fn a_homopolymer_of_six_or_more_is_a_period_one_locus_at_default() {
        let mut bases = filler(300);
        // 20 bp of poly-A, clean aperiodic flanks ≥ 30 bp either side (the default
        // flank is 30). At `--min-period 1` this is a period-1 locus.
        bases[100..120].copy_from_slice(&[b'A'; 20]);
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "homopolymer at default");

        let loci: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) => Some(l),
                _ => None,
            })
            .collect();
        assert_eq!(
            loci.len(),
            1,
            "the poly-A is a single period-1 locus at Default: {:#?}",
            kinds_of(&regions)
        );
        assert_eq!(loci[0].period(), 1, "period 1 — a homopolymer");
        assert_eq!(loci[0].motif().as_bytes(), b"A", "the homopolymer's motif");
        assert!(
            loci[0].start() >= 100 && loci[0].start() <= 105,
            "the locus is the poly-A run at ~101 (1-based), not a cascade artefact: {}",
            loci[0].start()
        );
    }

    /// A contig carrying exactly one repeat run, well clear of the flank
    /// requirement, with the bases **immediately either side forced** to one the
    /// run is not made of.
    ///
    /// The guard is the point. `filler` contains isolated `A`s, so dropping a
    /// five-base poly-A into it can land beside one and become a *six*-base run —
    /// and a floor test whose run is not the length it claims passes for the
    /// wrong reason. Forcing the neighbours makes the run's length exactly what
    /// the caller asked for.
    fn contig_with_one_guarded_run(run: &[u8], guard: u8) -> Vec<u8> {
        const OFFSET: usize = 200;
        assert!(
            !run.contains(&guard),
            "the guard must not be a base the run is made of, or it could extend it"
        );
        let mut bases = filler(OFFSET * 2 + run.len());
        bases[OFFSET - 1] = guard;
        bases[OFFSET..OFFSET + run.len()].copy_from_slice(run);
        bases[OFFSET + run.len()] = guard;
        bases
    }

    /// The `SsrSegment`s a default-config walk finds in `bases`.
    fn loci_at_default(bases: &[u8], case: &str) -> Vec<SsrSegment> {
        let regions = partition_resident("chr1", ContigId(0), bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), bases.len() as u64, case);
        regions
            .into_iter()
            .filter_map(|r| match r.kind {
                RegionKind::SsrSegment(l) => Some(l),
                _ => None,
            })
            .collect()
    }

    /// **The mononucleotide floor is exactly 8**: seven copies fall through as
    /// `Generic`, eight are a period-1 locus.
    ///
    /// `a_homopolymer_of_six_or_more_is_a_period_one_locus_at_default` uses a
    /// 20 bp run — it proves homopolymers are typed at all, but sits far from the
    /// line and so would stay green through any floor edit. This pins the line
    /// itself.
    ///
    /// **8 since 2026-08-10, measured, replacing an unmeasured 6.** The per-read
    /// slippage rate fitted over 181 tomato libraries first reaches 5% at 9 to 13
    /// repeats depending on the library; 8 takes one repeat of margin, since the
    /// most-stuttering library reaches 3.6% there against 1.0% at 7
    /// ([`MinCopies::default`](segment_criteria::MinCopies)).
    #[test]
    fn the_mononucleotide_copy_floor_is_exactly_eight() {
        let below = loci_at_default(&contig_with_one_guarded_run(&[b'A'; 7], b'C'), "poly-A x7");
        assert!(
            below.is_empty(),
            "7 copies is under the mono floor of 8 — the run stays Generic, got {below:?}"
        );

        let at = loci_at_default(&contig_with_one_guarded_run(&[b'A'; 8], b'C'), "poly-A x8");
        assert_eq!(at.len(), 1, "8 copies clears the mono floor: {at:?}");
        assert_eq!(at[0].period(), 1, "a homopolymer is period 1");
        assert_eq!(at[0].motif().as_bytes(), b"A");
        assert_eq!(at[0].tract_len(), 8, "the whole run, and only the run");
    }

    /// **The dinucleotide floor is exactly 6**: five copies fall through, six are a
    /// locus.
    ///
    /// **6 since 2026-08-10, and it is the one floor two independent measurements
    /// agree on.** The archive survey of 2,457 libraries puts the model-fit
    /// crossing at 6 — below it, most of what differs from the reference differs by
    /// something other than a whole number of copies — and the fitted slippage rate
    /// over 181 libraries first reaches 5% at 6 to 9 repeats. Both land on 6 for the
    /// most-stuttering library, which is the one a floor has to track.
    #[test]
    fn the_dinucleotide_copy_floor_is_exactly_six() {
        let below = loci_at_default(&contig_with_one_guarded_run(b"ATATATATAT", b'C'), "(AT)x5");
        assert!(
            below.is_empty(),
            "5 copies is under the di floor of 6 — the run stays Generic, got {below:?}"
        );

        let at = loci_at_default(
            &contig_with_one_guarded_run(b"ATATATATATAT", b'C'),
            "(AT)x6",
        );
        assert_eq!(at.len(), 1, "6 copies clears the di floor: {at:?}");
        assert_eq!(at[0].period(), 2);
        assert_eq!(at[0].motif().as_bytes(), b"AT");
        assert_eq!(at[0].tract_len(), 12, "six copies of a 2 bp motif");
    }

    /// **A rejected repeat is generic territory, not a hole** (spec §2.2). A
    /// low-copy tract classification turns down must still be *covered* — completeness
    /// is what the invariant is for.
    #[test]
    fn a_rejected_repeat_leaves_no_hole() {
        // (AT)*3 — below the period-2 copy floor (4 at the short-read Default).
        let mut bases = vec![b'C'; 60];
        bases.extend_from_slice(b"ATATAT");
        bases.extend(std::iter::repeat_n(b'G', 60));
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "rejected repeat");
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
            "below the copy floor: not a locus"
        );
        assert_eq!(
            regions.len(),
            1,
            "and the whole contig is Generic — no hole"
        );
    }

    /// A tract at position 1 has no left flank, so it is not genotypeable and lands
    /// in `Generic` (spec §8's edge list). The partition still tiles from base 1.
    #[test]
    fn a_tract_at_position_one_is_generic_and_the_partition_still_starts_at_one() {
        let mut bases = b"ATATATATATATATAT".to_vec();
        bases.extend(std::iter::repeat_n(b'C', 80));
        let len = bases.len() as u64;

        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "tract at base 1");
        assert_eq!(regions[0].region.start, Position(1));
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
            "no left flank to anchor against: not a locus"
        );
    }

    #[test]
    fn an_empty_contig_has_no_regions() {
        assert!(
            partition_resident("chr1", ContigId(0), b"", &TypedRegionConfig::default()).is_empty()
        );
    }

    // ---- D2: bundle detection --------------------------------------------

    /// Build a contig with `(AT)*10` tracts at the given 0-based offsets, aperiodic
    /// filler between and around them. Each tract is 20 bp.
    fn contig_with_tracts_at(offsets: &[usize], total: usize) -> Vec<u8> {
        let mut bases = filler(total);
        for &off in offsets {
            bases[off..off + 20].copy_from_slice(b"ATATATATATATATATATAT");
        }
        bases
    }

    /// **Two tracts 10 bp apart are ONE `SsrBundle` carrying both** — not two
    /// `Generic` regions, and not a hole. Neither can have a clean flank, so
    /// neither is a locus; but they are real repeats, and production would simply
    /// have deleted them (spec §2.4).
    #[test]
    fn two_close_tracts_become_one_bundle_carrying_both() {
        // Tracts at [60,80) and [90,110): a 10 bp gap, well inside bundle_threshold = 50.
        let bases = contig_with_tracts_at(&[60, 90], 240);
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "two close tracts");

        let bundles: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts)),
                _ => None,
            })
            .collect();
        assert_eq!(bundles.len(), 1, "ONE bundle, not two Generic regions");
        assert_eq!(bundles[0].1.len(), 2, "carrying BOTH tracts");
        // **The hull covers both tracts** — asserted as bounds, not as an exact
        // coordinate. The detector decides where a repeat starts and stops, and a
        // ±1–2 bp boundary/phase wobble is expected of it (`scanner_parity`
        // documents the same thing against trf-mod); pinning the exact edge would
        // be testing the detector's phase, not the bundle's hull.
        let (hull, tracts) = (bundles[0].0, bundles[0].1);
        assert!(hull.start <= Position(61), "hull reaches the first tract");
        assert!(hull.end >= Position(110), "hull reaches the last tract");
        assert_eq!(
            hull.start.get(),
            tracts.iter().map(|t| t.start).min().unwrap() + 1,
            "the hull IS the tracts' span, 1-based"
        );
        assert_eq!(hull.end.get(), tracts.iter().map(|t| t.end).max().unwrap());
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrSegment(_))),
            "neither tract has a clean flank, so neither is a locus"
        );
    }

    /// **Three tracts chained 30 bp apart are ONE bundle of three** — emergent
    /// transitivity. There is no separate transitive rule to implement: membership
    /// is the local flank test, and the cluster falls out of it (spec §2.4).
    ///
    /// A–B–C where A and C are 70 bp apart — further than `bundle_threshold` — still chain,
    /// because B is close to both.
    #[test]
    fn three_chained_tracts_become_one_bundle_of_three() {
        // [60,80), [110,130), [160,180): each gap 30 bp; A→C is 80 bp apart.
        let bases = contig_with_tracts_at(&[60, 110, 160], 300);
        let len = bases.len() as u64;
        let regions = partition_resident("chr1", ContigId(0), &bases, &catalog_config());
        assert_partitions(&regions, ContigId(0), len, "three chained tracts");

        let bundles: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts)),
                _ => None,
            })
            .collect();
        assert_eq!(bundles.len(), 1, "the chain is ONE bundle, not two");
        assert_eq!(
            bundles[0].1.len(),
            3,
            "all three, though the outer pair is further apart than bundle_threshold — \
             transitivity is emergent, not a rule"
        );
        // Bounds, not exact edges — the detector's ±1–2 bp phase wobble is its
        // business, not this test's.
        assert!(bundles[0].0.start <= Position(61));
        assert!(
            bundles[0].0.end >= Position(180),
            "the hull spans the whole chain, outermost tract to outermost tract"
        );
    }

    /// **Bundles do not spread.** A tract far enough from a cluster is classified
    /// normally — the flank test is local, so a bundle does not contaminate its
    /// neighbourhood.
    #[test]
    fn a_tract_clear_of_a_bundle_is_still_admitted() {
        // A close pair at [60,80)+[90,110), and a lone tract at [400,420) — far
        // from everything, with room for its 50 bp flanks.
        let bases = contig_with_tracts_at(&[60, 90, 400], 600);
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        assert_partitions(&regions, ContigId(0), len, "bundle + lone tract");

        assert_eq!(
            regions
                .iter()
                .filter(|r| matches!(r.kind, RegionKind::SsrBundle { .. }))
                .count(),
            1,
            "the close pair bundles"
        );
        let loci: Vec<_> = regions
            .iter()
            .filter_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) => Some(l),
                _ => None,
            })
            .collect();
        assert_eq!(
            loci.len(),
            1,
            "the lone tract is a locus — bundles don't spread"
        );
        assert_eq!(loci[0].start(), 401);
    }

    /// **A repeat rejected for any reason OTHER than the flank test is `Generic`,
    /// not a bundle** (spec §2.2). Only closeness makes a bundle; being low-copy
    /// makes you ordinary.
    #[test]
    fn a_low_copy_repeat_is_generic_not_a_bundle() {
        let mut bases = filler(240);
        // (AT)*3 — below the period-2 copy floor (4 at the short-read Default), and isolated.
        bases[100..106].copy_from_slice(b"ATATAT");
        let len = bases.len() as u64;
        let regions =
            partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());

        assert_partitions(&regions, ContigId(0), len, "low-copy repeat");
        assert!(
            !regions
                .iter()
                .any(|r| matches!(r.kind, RegionKind::SsrBundle { .. })),
            "a copy-floor rejection is Generic territory, NOT a bundle"
        );
        assert_eq!(regions.len(), 1, "and no hole: the contig is one Generic");
    }

    /// The clustering must reproduce exactly the split that produced it — a
    /// singleton "bundle" would mean the two disagree.
    #[test]
    fn bundle_clusters_regroups_exactly_what_the_split_set_aside() {
        let iv = |start, end| RepeatInterval {
            start,
            end,
            period: 2,
            score: 100,
        };
        // Two clusters: {100,130 / 150,180} and {5000,5030 / 5040,5070}.
        let bundled = [iv(100, 130), iv(150, 180), iv(5000, 5030), iv(5040, 5070)];
        let clusters = segment_criteria::bundle_clusters(&bundled, 50);
        assert_eq!(clusters.len(), 2, "two clusters, not one and not four");
        assert_eq!(clusters[0].len(), 2);
        assert_eq!(clusters[1].len(), 2);
        assert_eq!(clusters[0][0].start, 100);
        assert_eq!(clusters[1][0].start, 5000);
    }

    /// **`.cat` parity — the oracle, and D1's real bar** (spec §8.1).
    ///
    /// The walk at the catalog's settings must reproduce the catalog: every golden
    /// locus is present, **or absent *and* inside a satellite run**. A strict
    /// subset, and that shape is *earned* by spec §2.4's ordering — the cap applies
    /// to the *cleaned* coverage, after classification, so the difference can only go
    /// one way. Capping raw coverage would make it bidirectional and untestable.
    ///
    /// **The oracle is the committed trf-mod-built golden catalog** — a different
    /// detector, a different code path, nothing ng touched. Its detector difference
    /// is already characterised (`scanner_parity`: 16/16, 15 exact, one ±1–2 bp
    /// boundary/phase wobble, one genuine scanner-only locus trf-mod's significance
    /// model rejected), so it is a yardstick rather than a confound — hence overlap
    /// matching, inherited from `scanner_parity` for exactly that reason.
    ///
    /// What this proves: the plumbing — the scan, the pre-filter, the classification
    /// call, the coordinate arithmetic. What it does **not** prove: that the
    /// catalog's settings are *right* (spec §5). A fixed-config regression test,
    /// pinned to those settings explicitly rather than to whatever `Default` is, or
    /// it starts failing the first time someone moves a floor and reads a *result*
    /// as a bug.
    #[test]
    fn the_resident_partition_reproduces_the_golden_catalog() {
        use crate::ssr::catalog::io::CatalogReader;
        use std::fs::File;
        use std::io::BufReader;
        use std::path::Path;

        let fixture = |name: &str| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests")
                .join("data")
                .join("tandem_repeat")
                .join(name)
        };

        // The golden catalog and the settings it was built at — read, never
        // written (production is frozen).
        let mut golden_reader =
            CatalogReader::new(File::open(fixture("golden.ssr_catalog.bed.gz")).unwrap()).unwrap();
        let cat_params = golden_reader.header().params.clone();
        let golden = golden_reader.read_all().unwrap();
        assert!(!golden.is_empty(), "the golden catalog must have loci");

        // ng's walk at the SAME settings, pinned explicitly. As of spec §2.3 ng's
        // `Default` carries the short-read floors, not the catalog's, so the period
        // scope, copy floors, and satellite cap are pinned to the catalog's build
        // settings here (§8.1); `min_purity`/`min_score`/`bundle_threshold` come from the
        // `.cat` header.
        let config = TypedRegionConfig {
            criteria: SsrSegmentCriteria {
                min_purity: cat_params.min_purity,
                min_score: cat_params.min_score,
                bundle_threshold: u64::from(cat_params.bundle_threshold),
                periods: crate::ng::tandem_repeat::PeriodRange::new(2, 6)
                    .expect("2..=6 is a valid period range"),
                min_copies: segment_criteria::MinCopies::new([10, 5, 4, 3, 3, 3], 3),
            },
            max_str_len: Bp(1000),
            ..TypedRegionConfig::default()
        };

        let file = File::open(fixture("synthetic_ref.fa")).unwrap();
        let mut reader = noodles_fasta::io::Reader::new(BufReader::new(file));
        let mut ours: Vec<(String, u64, u64)> = Vec::new();
        let mut satellites: Vec<(String, u64, u64)> = Vec::new();
        for (idx, result) in reader.records().enumerate() {
            let rec = result.unwrap();
            let name = String::from_utf8_lossy(rec.name()).into_owned();
            let bases = rec.sequence().as_ref();
            let regions = partition_resident(&name, ContigId(idx as u32), bases, &config);

            // The partition must hold on REAL sequence, not only the crafted
            // fixtures above.
            assert_partitions(
                &regions,
                ContigId(idx as u32),
                bases.len() as u64,
                &format!("golden contig {name}"),
            );

            for r in &regions {
                match &r.kind {
                    RegionKind::SsrSegment(l) => ours.push((name.clone(), l.start(), l.end())),
                    RegionKind::Satellite => {
                        satellites.push((name.clone(), r.region.start.get(), r.region.end.get()))
                    }
                    _ => {}
                }
            }
        }
        assert!(!ours.is_empty(), "the scan must find loci");

        // Overlap match. Production's `Locus` is 0-based half-open, ng's is 1-based
        // inclusive, so production's `[s, e)` is ng's `[s + 1, e]` (spec §4).
        let overlaps =
            |a: &(String, u64, u64), b: &(String, u64, u64)| a.0 == b.0 && a.1 <= b.2 && b.1 <= a.2;
        let mut missed = Vec::new();
        for g in &golden {
            let g1 = (
                g.chrom().to_string(),
                u64::from(g.start()) + 1,
                u64::from(g.end()),
            );
            if !ours.iter().any(|o| overlaps(&g1, o)) {
                // Absent is legal ONLY inside a satellite run — the one expected
                // divergence (the file applies no cap; a reader does).
                if !satellites.iter().any(|s| overlaps(&g1, s)) {
                    missed.push(format!("{}:{}-{}", g1.0, g1.1, g1.2));
                }
            }
        }
        assert!(
            missed.is_empty(),
            "every golden locus must be present, or absent AND inside a satellite \
             run. At the catalog's settings a locus missing for any other reason is \
             a machinery bug. Missing: {missed:#?}"
        );
    }

    // ---- D3: the windowed walk -------------------------------------------

    fn windowing_fixture() -> Vec<u8> {
        let mut bases = filler(6000);
        let mut tract_at = |off: usize| {
            bases[off..off + 20].copy_from_slice(b"ATATATATATATATATATAT");
        };
        tract_at(990);
        tract_at(1990);
        tract_at(2020);
        tract_at(3500);
        for i in 0..600 {
            bases[4000 + i * 2..4002 + i * 2].copy_from_slice(b"AT");
        }
        bases
    }

    /// A **chain**: 20 bp `(AT)` tracts every 50 bp — a 30 bp gap between each, under the
    /// 50 bp flank — running from 1000 for `span_bp`. Every tract is close to the next, so
    /// the whole thing is one cluster however long it is.
    fn chained_cluster(span_bp: usize, total: usize) -> Vec<u8> {
        let mut bases = filler(total);
        let mut at = 1000;
        while at + 20 <= 1000 + span_bp {
            bases[at..at + 20].copy_from_slice(b"ATATATATATATATATATAT");
            at += 50;
        }
        bases
    }

    /// **A BED returns the same findings a whole-genome run does — the same *objects*, not
    /// merely the same kinds** (spec §2.5, owner 2026-07-17). This is what the whole-contig
    /// scan set buys, and it is the strongest form of BED-invariance available.
    ///
    /// # The two fixtures that broke the old design
    ///
    /// E2 grew each requested span by `max_str_len` and scanned that. Both of these
    /// defeat a margin, and both are pinned here because they are exactly the shapes a
    /// margin cannot hold:
    ///
    /// - **a 3 kb array** — a satellite is *by definition* longer than `max_str_len`,
    ///   and the margin **is** `max_str_len`. Under the old design a BED reported it as
    ///   `1001–2300` where the truth is `1001–4001`: a `Satellite` cut to 1300 bp,
    ///   silently. This is the common case, not an exotic one.
    /// - **a 1.5 kb chain** of tracts 30 bp apart — clustering has no reach at all (spec
    ///   §2.6: *"A–B–C–D each 30 bp apart runs past any margin you choose"*). §2.5 hoped
    ///   such a chain would be *"dense-repeat territory heading for `Satellite`"*; it is
    ///   not — a satellite is over-cap **coverage**, and coverage runs merge only where
    ///   they **abut**, so every run in this chain is 20 bp. It stays a bundle, and the old
    ///   design cut it to 24 members of 30.
    #[test]
    fn a_bed_returns_the_same_findings_the_whole_genome_run_does() {
        // The satellite-chain fixtures below are tuned to the catalog's 50 bp
        // bundle radius and 1 kb cap, so this pins them explicitly (spec §2.3
        // moved `Default` off those values). The BED-invariance property under
        // test is config-independent; the fixture geometry is not.
        let config = catalog_config();

        // --- the 3 kb array: a satellite three times the old margin ---
        let bases = contig_with_a_long_array();
        let whole = partition_resident("chr1", ContigId(0), &bases, &config);
        let truth = whole
            .iter()
            .find(|r| matches!(r.kind, RegionKind::Satellite))
            .expect("the fixture has one array");
        assert!(
            truth.region.len() > 3000,
            "the array must be far longer than max_str_len: {:?}",
            truth.region
        );

        // Ask about 200 bp of its left end. The satellite comes back WHOLE — all 3 kb,
        // past both requested edges, identical to the whole-genome region.
        let got = walk_bed_with(&bases, &[(1100, 1300)], catalog_config());
        let sat = got
            .iter()
            .find(|r| matches!(r.kind, RegionKind::Satellite))
            .expect("still a satellite");
        assert_eq!(
            sat.region, truth.region,
            "a BED must not shorten an array: the extent IS the claim, and a margin cannot \
             hold a feature defined as being longer than it"
        );

        // --- the 1.5 kb chain: a bundle no margin could hold ---
        let bases = chained_cluster(1500, 6000);
        let whole = partition_resident("chr1", ContigId(0), &bases, &config);
        let (truth_hull, truth_tracts) = whole
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts.clone())),
                _ => None,
            })
            .expect("the chain is one bundle");
        assert!(
            truth_hull.len() > 1400 && truth_tracts.len() >= 25,
            "the chain must run well past max_str_len: {truth_hull:?}"
        );
        assert!(
            !whole
                .iter()
                .any(|r| matches!(r.kind, RegionKind::Satellite)),
            "and it is NOT a satellite — §2.5 hoped it would be, but coverage runs merge \
             only where they abut and every run here is 20 bp: {whole:#?}"
        );

        let got = walk_bed_with(&bases, &[(900, 1200)], catalog_config());
        let (hull, tracts) = got
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrBundle { tracts } => Some((r.region, tracts.clone())),
                _ => None,
            })
            .expect("still a bundle");
        assert_eq!(hull, truth_hull, "the same hull");
        assert_eq!(
            tracts, truth_tracts,
            "and the same member tracts — all 30, not the 24 inside a grown span"
        );
    }

    /// A `(CAG)*8` microsatellite `gap` bp from a 1.3 kb `(AT)` array, on either side of
    /// it. Different motifs, so the detector cannot join the two into one tract.
    fn micro_near_satellite(micro_left: bool, gap: usize) -> Vec<u8> {
        let mut bases = filler(6000);
        let (micro_at, array_at) = if micro_left {
            (1000, 1000 + 24 + gap)
        } else {
            (1000 + 1300 + gap, 1000)
        };
        // (CAG)*8 — a different motif from the array's, so the detector cannot join them.
        for i in 0..8 {
            bases[micro_at + i * 3..micro_at + i * 3 + 3].copy_from_slice(b"CAG");
        }
        for i in 0..650 {
            bases[array_at + i * 2..array_at + i * 2 + 2].copy_from_slice(b"AT");
        }
        bases
    }

    /// **A microsatellite too close to a satellite is swallowed by it, from either
    /// side** (spec §2.4a — the owner's rule, 2026-07-16).
    ///
    /// The situation is physical and symmetric: a `(CAG)*8` tract close to a 1.3 kb
    /// array cannot be genotyped, because the flank on that side *is* array. So the
    /// answer must not depend on which side it sits — and before this rule it did, in
    /// two different wrong ways (`resolve_features`' docs).
    ///
    /// The far arm is the **control**, and it is what makes the near arm mean
    /// something: 200 bp away the same two features, built by the same code, give a
    /// clean `Generic / locus / Generic / Satellite` — so absorption up close is the
    /// *rule* firing, not classification quietly rejecting a tract that was never
    /// viable.
    ///
    /// **Both gaps are derived from [`DEFAULT_BUNDLE_THRESHOLD`] rather than written as
    /// numbers.** They used to be a literal 20 bp against a radius of 30; when the
    /// radius moved to 20 the near arm landed exactly on the boundary — where the rule
    /// is strict and so does *not* absorb — and the test failed for a reason that had
    /// nothing to do with what it checks. A fixture that pins a constant it does not
    /// own tests the constant.
    #[test]
    fn a_microsatellite_beside_a_satellite_is_absorbed_into_it() {
        let config = TypedRegionConfig::default();
        // Comfortably inside the radius, and comfortably outside it.
        let near = (DEFAULT_BUNDLE_THRESHOLD / 2) as usize;
        let far = (DEFAULT_BUNDLE_THRESHOLD * 10) as usize;
        let kinds = |regions: &[TypedRegion]| -> Vec<&'static str> {
            regions
                .iter()
                .map(|r| match &r.kind {
                    RegionKind::SsrSegment(_) => "locus",
                    RegionKind::SsrBundle { .. } => "bundle",
                    RegionKind::Generic => "generic",
                    RegionKind::Satellite => "satellite",
                })
                .collect()
        };

        for micro_left in [true, false] {
            let side = if micro_left { "left" } else { "right" };

            // --- absorbed: a gap inside the bundle radius ---
            let bases = micro_near_satellite(micro_left, near);
            let regions = partition_resident("chr1", ContigId(0), &bases, &config);
            assert_partitions(
                &regions,
                ContigId(0),
                bases.len() as u64,
                &format!("micro {side} of a satellite, {near} bp"),
            );
            assert_eq!(
                kinds(&regions),
                vec!["generic", "satellite", "generic"],
                "micro {side}: the satellite swallows it — no bundle, no locus, and \
                 (crucially) no second region over the same bases"
            );

            // The satellite **grew to cover the microsatellite**: it is not merely that
            // the tract went missing. The absorbed span reaches past the array's own
            // 1.3 kb, gap included.
            let satellite = regions
                .iter()
                .find(|r| matches!(r.kind, RegionKind::Satellite))
                .expect("one satellite");
            assert!(
                satellite.region.len() as usize >= 1300 + near + 24,
                "micro {side}: the satellite must EXPAND over the gap and the tract — \
                 array (1300) + gap ({near}) + tract (24); got {}",
                satellite.region.len()
            );

            // --- the control: far outside the radius, and the flank is clean ---
            let bases = micro_near_satellite(micro_left, far);
            let regions = partition_resident("chr1", ContigId(0), &bases, &config);
            assert_partitions(
                &regions,
                ContigId(0),
                bases.len() as u64,
                &format!("micro {side} of a satellite, {far} bp"),
            );
            assert!(
                kinds(&regions).contains(&"locus"),
                "micro {side}: {far} bp away the SAME tract is a locus — so the absorption \
                 above is the rule, not a tract that was never admissible"
            );
            assert!(kinds(&regions).contains(&"satellite"));
        }
    }

    /// **A wider radius absorbs from further away.** That is the whole content of the
    /// knob, and it is what nothing asserted until a `(CAG)*8` tract beside a 1.3 kb array
    /// started behaving differently when the default moved.
    ///
    /// Pinned because its absence cost two confusing failures. The absorption test above
    /// used a literal 20 bp gap, comfortable inside a radius of 30 and right at the
    /// transition once the radius became 20. Then a first version of *this* test asserted
    /// "absorbed at exactly the radius, a locus one base further" — true at a radius of 20
    /// and false at 15. The behaviour was right every time; the assertions were guesses at
    /// a relationship nobody had measured.
    ///
    /// **So it asserts the monotonicity and not an offset.** How far a tract must sit from
    /// the array before it survives is not the radius plus a constant: classification trims
    /// a tract's edges, so the gap this fixture lays down and the gap the rule compares are
    /// not the same number, and the difference moves with the radius. What is guaranteed,
    /// and what a reader needs, is the direction — and it is checked across four radii
    /// rather than at the default alone, so moving the default cannot break it again.
    #[test]
    fn absorption_tracks_the_bundle_radius() {
        let first_locus_gap = |radius: u64| -> Option<usize> {
            let config = TypedRegionConfig {
                criteria: SsrSegmentCriteria {
                    bundle_threshold: radius,
                    ..SsrSegmentCriteria::default()
                },
                ..TypedRegionConfig::default()
            };
            (1..=radius as usize * 3).find(|&gap| {
                let bases = micro_near_satellite(true, gap);
                partition_resident("chr1", ContigId(0), &bases, &config)
                    .iter()
                    .any(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
            })
        };
        let mut previous = 0;
        for radius in [10u64, 15, 20, 30] {
            let gap = first_locus_gap(radius)
                .unwrap_or_else(|| panic!("radius {radius}: no gap leaves a locus at all"));
            assert!(
                gap > previous,
                "radius {radius}: the tract should survive at a wider gap than the previous \
                 radius did — got {gap} against {previous}"
            );
            previous = gap;
        }
    }

    /// **The cap applies to the CLEANED coverage, not the raw** (spec §2.4).
    ///
    /// The noise that discriminates the two sets is now the **copy floor's**: since the
    /// scanner emits only primitive periods, its raw output no longer carries the flood
    /// of aliases D1 relied on (raw and cleaned coverage coincide on the plain windowing
    /// fixture — 5 runs each). What still differs is a *low-copy* interval: the scanner
    /// emits down to 2 copies, and the copy floor drops it. Here a 2-copy `(CG)` speck
    /// abuts the 1.2 kb array's right edge — in the RAW coverage it extends the array's
    /// run 4 bp past the cap; in the CLEANED coverage it is gone. Cap the wrong set and
    /// the satellite's edge moves.
    ///
    /// Written against the two candidate computations rather than a coordinate literal:
    /// it must fail when the cap moves to the raw set, and **not** when the detector's
    /// phase shifts by a base, which is the detector's business (`scanner_parity`).
    ///
    /// **What this does not do.** The difference it catches is a boundary, not the
    /// failure the spec's argument is about — noise unioning past 1 kb, declaring a
    /// satellite of its own, and swallowing the loci underneath. No fixture produces
    /// that. This pins the choice; it does not demonstrate the stakes.
    ///
    /// **The catalog's reader shares this decision rather than repeating it** — it calls
    /// [`coverage_runs`] over its own cleaned set — so this test covers both paths, and
    /// the differential cannot see the choice move because it would move on both sides at
    /// once.
    #[test]
    fn the_satellite_cap_applies_to_the_cleaned_coverage_not_the_raw() {
        let mut bases = windowing_fixture();
        // A 2-copy (CG) speck abutting the 1.2 kb (AT) array's right edge (it ends at
        // 5200). The scanner emits it (min_copies 2); the copy floor (period 2 → 5)
        // drops it. So the raw coverage run runs 4 bp longer than the cleaned one.
        bases[5200..5204].copy_from_slice(b"CGCG");
        let config = TypedRegionConfig::default();
        let raw = find_tandem_repeats(&bases, config.criteria.periods, &config.scan);
        let cleaned = segment_criteria::prefilter(&raw, &config.criteria);

        let over_cap = |intervals: &[RepeatInterval]| -> Vec<(u64, u64)> {
            coverage_runs(intervals)
                .into_iter()
                .filter(|r| r.len() > config.max_str_len.get())
                .map(|r| (r.start, r.end))
                .collect()
        };
        let from_raw = over_cap(&raw);
        let from_cleaned = over_cap(&cleaned);
        assert_ne!(
            from_raw, from_cleaned,
            "this fixture must DISCRIMINATE the two, or the assertion below cannot fail \
             — which is exactly how the claim went untested through D1 and D2"
        );

        let satellites: Vec<(u64, u64)> = partition_resident("chr1", ContigId(0), &bases, &config)
            .iter()
            .filter(|r| matches!(r.kind, RegionKind::Satellite))
            .map(|r| (r.region.start.get(), r.region.end.get()))
            .collect();
        assert_eq!(
            satellites, from_cleaned,
            "the satellites are the CLEANED coverage's over-cap runs: capping the \
             raw set would let detector noise decide where an array begins — and, with 1 \
             kb of it in a row, that an array exists at all (spec §2.4)"
        );
    }

    // ---- The tally -------------------------------------------------------

    /// The **tally describes the regions that came out** (spec §3.1).
    ///
    /// Checked against the regions actually yielded rather than against literals — the
    /// counts are a claim *about the output*, so anything else would be two independent
    /// guesses at the fixture.
    #[test]
    fn counts_tally_the_regions_yielded() {
        let bases = windowing_fixture();
        let (yielded, counts) = partition_resident_in(
            "chr1",
            ContigId(0),
            &bases,
            &TypedRegionConfig::default(),
            &whole_contig(bases.len() as u64),
        );
        let counts = &counts;
        let kind_count =
            |f: fn(&RegionKind) -> bool| yielded.iter().filter(|r| f(&r.kind)).count() as u64;
        assert_eq!(counts.spans, 1);
        assert_eq!(
            counts.ssr_loci,
            kind_count(|k| matches!(k, RegionKind::SsrSegment(_)))
        );
        assert_eq!(
            counts.ssr_bundles,
            kind_count(|k| matches!(k, RegionKind::SsrBundle { .. }))
        );
        assert_eq!(
            counts.generic,
            kind_count(|k| matches!(k, RegionKind::Generic))
        );
        assert_eq!(
            counts.satellites,
            kind_count(|k| matches!(k, RegionKind::Satellite))
        );
        assert!(counts.satellites > 0 && counts.ssr_loci > 0 && counts.ssr_bundles > 0);

        let bp_of = |f: fn(&RegionKind) -> bool| -> u64 {
            yielded
                .iter()
                .filter(|r| f(&r.kind))
                .map(|r| r.region.len())
                .sum()
        };
        assert_eq!(
            counts.satellite_bp,
            bp_of(|k| matches!(k, RegionKind::Satellite))
        );
        assert_eq!(
            counts.ssr_bundle_bp,
            bp_of(|k| matches!(k, RegionKind::SsrBundle { .. }))
        );

        // **The gap `repeat_bp_with_no_locus` measures**, pinned exactly rather than
        // bounded: it is the contig's cleaned repeat coverage minus the bases that came
        // out as loci. Computed here from the resident path — an independent route to the
        // same number, so this fails if either the accumulate or the subtract is wrong.
        // (A `>=` bound was the first version, and it left the subtraction untested:
        // dropping it makes the number bigger, and bigger still satisfies `>=`.)
        let config = TypedRegionConfig::default();
        let raw = find_tandem_repeats(&bases, config.criteria.periods, &config.scan);
        let cleaned = segment_criteria::prefilter(&raw, &config.criteria);
        let coverage_bp: u64 = coverage_runs(&cleaned).iter().map(|r| r.len()).sum();
        let locus_bp = bp_of(|k| matches!(k, RegionKind::SsrSegment(_)));
        assert!(
            coverage_bp > 0 && locus_bp > 0,
            "the fixture must have both"
        );
        assert_eq!(
            counts.repeat_bp_with_no_locus,
            coverage_bp - locus_bp,
            "repeat coverage ({coverage_bp} bp) that yielded no locus ({locus_bp} bp of \
             it did)"
        );
        assert!(
            counts.repeat_bp_with_no_locus >= counts.ssr_bundle_bp,
            "and the bundled tracts are part of it"
        );
    }

    /// One contig, end to end, as the region list [`partition_resident_in`] is asked with.
    fn whole_contig(len: u64) -> [GenomeRegion; 1] {
        [GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(len),
        }]
    }

    /// The finished tally over a whole contig, at an explicit config.
    fn counts_over_with(bases: &[u8], config: TypedRegionConfig) -> TypedRegionCounts {
        partition_resident_in(
            "chr1",
            ContigId(0),
            bases,
            &config,
            &whole_contig(bases.len() as u64),
        )
        .1
    }

    // ---- Asking for part of a contig -------------------------------------

    /// The regions of `bases` that fall inside the 1-based inclusive spans `want`, at an
    /// explicit config. `walk_bed` is this at `Default`.
    fn walk_bed_with(
        bases: &[u8],
        want: &[(u64, u64)],
        config: TypedRegionConfig,
    ) -> Vec<TypedRegion> {
        let requested: Vec<GenomeRegion> = want
            .iter()
            .map(|&(s, e)| GenomeRegion {
                contig: ContigId(0),
                start: Position(s),
                end: Position(e),
            })
            .collect();
        partition_resident_in("chr1", ContigId(0), bases, &config, &requested).0
    }

    fn walk_bed(bases: &[u8], want: &[(u64, u64)]) -> Vec<TypedRegion> {
        walk_bed_with(bases, want, TypedRegionConfig::default())
    }

    /// **BED-invariance — the acceptance test** (spec §2.5): *"whether a base is STR /
    /// bundle / generic / satellite must not depend on the BED. The BED chooses what you
    /// are shown, never what things are."*
    ///
    /// Every requested span is compared against the **whole-genome** run restricted to the
    /// same coordinates. The spans are chosen to cut through everything the fixture has —
    /// a BED edge 5 bp inside the bundle's cluster, one mid-satellite, one splitting a
    /// generic run — because an edge in the middle of nowhere tests nothing.
    ///
    /// This is what "scan wider than you emit" buys, and it is not free: the naive walk
    /// (scan only what was asked for) gets the *bundle* wrong, because the neighbour that
    /// makes it a bundle sits outside the span.
    #[test]
    fn a_bed_does_not_change_what_things_are() {
        let bases = windowing_fixture();
        let whole = partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());

        // Spans that cut through real structure: across the tract at 990; across the
        // bundle at 1990..2040; inside the satellite at 4000..5200; and a plain interior.
        for want in [
            (900u64, 1100u64),
            (1995, 2100),
            (4500, 4800),
            (3000, 3600),
            (1, 6000),
        ] {
            let got = walk_bed(&bases, &[want]);

            // What the whole-genome run says about those same bases.
            for r in &got {
                for pos in [r.region.start.get(), r.region.end.get()] {
                    // Only ask about bases the user actually requested: an object
                    // emitted whole reaches past the edge, which is the point.
                    if pos < want.0 || pos > want.1 {
                        continue;
                    }
                    let truth = whole
                        .iter()
                        .find(|w| w.region.contains(Position(pos)))
                        .unwrap_or_else(|| panic!("the whole-genome run covers {pos}"));
                    assert_eq!(
                        std::mem::discriminant(&r.kind),
                        std::mem::discriminant(&truth.kind),
                        "BED {want:?}: base {pos} is {:?} \
                             with the BED and {:?} without it — the BED changed what a \
                             base IS (spec §2.5)",
                        r.kind,
                        truth.kind
                    );
                }
            }

            // And the loci are the same objects, not merely the same kind.
            let loci_in = |rs: &[TypedRegion]| -> Vec<(u64, u64)> {
                rs.iter()
                    .filter(|r| matches!(r.kind, RegionKind::SsrSegment(_)))
                    .map(|r| (r.region.start.get(), r.region.end.get()))
                    .filter(|(s, _)| *s >= want.0 && *s <= want.1)
                    .collect()
            };
            assert_eq!(
                loci_in(&got),
                loci_in(&whole),
                "BED {want:?}: the loci starting inside the \
                     span must be exactly the whole-genome run's"
            );
        }
    }

    /// The requested span is **tiled exactly** — the partition invariant, restated where
    /// it belongs once a BED is involved: over the *effective* region, which is what the
    /// user asked for grown to hold a straddling object whole (spec §2.5).
    #[test]
    fn a_bed_span_is_tiled_by_what_comes_back() {
        let bases = windowing_fixture();
        for want in [(900u64, 1100u64), (1995, 2100), (3000, 3600)] {
            let got = walk_bed(&bases, &[want]);
            assert!(!got.is_empty());

            // Contiguous and non-overlapping, and covering every requested base. The
            // effective region may start before / end after the request (an object
            // straddling the edge), so the *bounds* are asserted as reaching at least the
            // request, not equalling it.
            assert!(got[0].region.start.get() <= want.0);
            assert!(got.last().unwrap().region.end.get() >= want.1);
            for pair in got.windows(2) {
                assert_eq!(
                    pair[1].region.start.get(),
                    pair[0].region.end.get() + 1,
                    "gap or overlap inside a BED span: {want:?}"
                );
                assert_ne!(
                    std::mem::discriminant(&pair[0].kind),
                    std::mem::discriminant(&pair[1].kind),
                    "two consecutive regions share a kind — MAXIMALITY (spec §2.3)"
                );
            }
        }
    }

    /// **`Generic` is clipped to the user's edge; every *finding* straddling it comes back
    /// whole** — locus, bundle **and satellite** (spec §2.5, owner 2026-07-17;
    /// [`clips_at_a_bed_edge`]).
    #[test]
    fn generic_clips_at_the_edge_and_findings_come_back_whole() {
        let bases = windowing_fixture();

        // A span whose edges fall in plain generic sequence: nothing straddles, so the
        // partition starts and ends exactly where the user asked.
        let got = walk_bed(&bases, &[(3000, 3600)]);
        assert_eq!(
            got[0].region.start,
            Position(3000),
            "Generic is clipped to the user's edge, not grown to the scan span"
        );
        assert!(matches!(got[0].kind, RegionKind::Generic));
        assert_eq!(got.last().unwrap().region.end, Position(3600));

        // A span wholly INSIDE the 1.2 kb array: the satellite comes back whole, reaching
        // past both requested edges.
        //
        // **Clipping it would produce a `Satellite` region of 301 bp** — a span that
        // contradicts the `max_str_len` (1 kb) test that produced the label. The extent
        // is the claim: "an array too long to be a microsatellite". That is what E2 got
        // wrong by reasoning from the type (`Satellite` carries no payload, so nothing
        // could be left misdescribed) instead of from the meaning.
        let got = walk_bed(&bases, &[(4500, 4800)]);
        assert_eq!(got.len(), 1, "the whole span is inside the array: {got:#?}");
        assert!(matches!(got[0].kind, RegionKind::Satellite));
        assert!(
            got[0].region.start.get() < 4500 && got[0].region.end.get() > 4800,
            "the Satellite is emitted WHOLE, past both edges: {:?}",
            got[0].region
        );
        assert!(
            got[0].region.len() > 1000,
            "and its span is over the cap that made it a satellite — which a clipped one \
             ({} bp of request) could not be",
            4800 - 4500 + 1
        );
        // The same region the whole-genome run reports, not a version of it.
        let whole = partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        let truth = whole
            .iter()
            .find(|r| matches!(r.kind, RegionKind::Satellite))
            .expect("the fixture has one");
        assert_eq!(got[0].region, truth.region, "the SAME satellite");

        // A span cutting into the bundle's cluster: the bundle carries its member tracts,
        // so clipping it would leave them outside their own region. It comes back whole,
        // and the effective region grows to hold it.
        let got = walk_bed(&bases, &[(1995, 2100)]);
        let bundle = got
            .iter()
            .find(|r| matches!(r.kind, RegionKind::SsrBundle { .. }))
            .expect("the cluster at 1990..2040 straddles this edge");
        assert!(
            bundle.region.start.get() < 1995,
            "the bundle is emitted WHOLE, reaching back past the requested edge to {} \
             (the effective region grows to hold it)",
            bundle.region.start.get()
        );
        let RegionKind::SsrBundle { tracts } = &bundle.kind else {
            unreachable!()
        };
        assert!(
            tracts
                .iter()
                .all(|t| t.start + 1 >= bundle.region.start.get()
                    && t.end <= bundle.region.end.get()),
            "every member tract is inside the hull — which is what clipping would break"
        );
    }

    /// A locus straddling the edge likewise comes back **whole**, with its bases intact —
    /// half a locus is not a locus.
    #[test]
    fn a_locus_straddling_the_bed_edge_is_emitted_whole() {
        let bases = windowing_fixture();
        // The tract at 990..1010 (1-based 991..1010); ask for a span ending inside it.
        let got = walk_bed(&bases, &[(900, 1000)]);
        let locus = got
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) => Some(l),
                _ => None,
            })
            .expect("the tract at ~991 straddles the requested edge at 1000");
        assert!(
            locus.end() > 1000,
            "the locus reaches past the user's edge ({}), whole",
            locus.end()
        );
        // The same object the whole-genome run builds, bases and all.
        let whole = partition_resident("chr1", ContigId(0), &bases, &TypedRegionConfig::default());
        let truth = whole
            .iter()
            .find_map(|r| match &r.kind {
                RegionKind::SsrSegment(l) if l.start() == locus.start() => Some(l),
                _ => None,
            })
            .expect("the same locus exists in the whole-genome run");
        assert_eq!(locus, truth, "and it is the SAME locus, not a clipped one");
    }

    /// **Two requested spans on one contig share its classification**: the contig is
    /// classified once, and the ground between them — which the caller did not ask for —
    /// must not come back (spec §2.5).
    ///
    /// This is the case that makes `Generic` clip against *each* requested span rather than
    /// against the scan span: one generic run covers both, and it has to come back as two
    /// regions with the gap dropped. Since the scan set became whole contigs it is no
    /// longer an edge case at all — **every** pair of spans on a contig is this case, which
    /// is a good reason for the rule to be the general one.
    #[test]
    fn two_spans_sharing_a_scan_span_do_not_leak_the_gap_between_them() {
        let bases = windowing_fixture();
        let got = walk_bed(&bases, &[(3000, 3200), (3400, 3600)]);

        assert!(
            got.iter()
                .all(|r| (3000..=3200).contains(&r.region.start.get())
                    || (3400..=3600).contains(&r.region.start.get())),
            "nothing may start in the gap the user did not ask for: {got:#?}"
        );
        let covered: u64 = got.iter().map(|r| r.region.len()).sum();
        assert_eq!(
            covered,
            201 + 201,
            "exactly the two requested spans come back — the gap between them is scanned \
             (they share a scan span) and not emitted"
        );
    }

    /// Every BED failure is `RegionSet`'s to reject **up front**, so a consumer holding
    /// a [`GenomeRegions`] has nothing left to validate (spec §8.2).
    #[test]
    fn a_bad_bed_is_rejected_before_any_region_is_typed() {
        use std::io::Write;
        std::fs::create_dir_all("tmp").unwrap();
        let dir = tempfile::tempdir_in("tmp").unwrap();
        let bed = dir.path().join("bad.bed");
        {
            let mut f = std::fs::File::create(&bed).unwrap();
            writeln!(f, "nosuchcontig\t0\t10").unwrap();
        }
        assert!(
            GenomeRegions::from_bed_path(&bed, CONTIGS).is_err(),
            "an unknown contig is caught at construction, not mid-walk"
        );
    }
}
