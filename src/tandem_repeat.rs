//! A clean-room, use-agnostic short-period **tandem-repeat scanner** — a sequence
//! primitive shared by the STR catalog and (later) the ng snp/str caller. Design:
//! spec `doc/devel/ng/spec/ssr_repeat_scanner.md`, arch
//! `doc/devel/ng/arch/ssr_repeat_scanner.md`.
//!
//! It exposes **three interfaces** over one detection core (a lag-`p` self-comparison
//! plus a Ruzzo–Tompa maximal-scoring-segment pass):
//!
//! - a low-level **interval finder** — `find_tandem_repeats(seq, PeriodRange,
//!   &ScanParams) -> Vec<RepeatInterval>` — raw, possibly-overlapping intervals, for
//!   consumers that resolve overlaps themselves (the STR catalog); and
//! - a mid-level **windowed scan** — `scan_windowed`, which yields one
//!   `ScannedWindow` at a time (the core + margin geometry, every detection in the
//!   fetched slice, and the two window rules — coverage clipped to cores, intervals
//!   attributed by start — as methods the consumer applies to *its own* interval set)
//!   and takes **no routing decision**; for a consumer that must inject its own policy
//!   between detection and any merge — the typed-region generator, which is why this
//!   layer is public (`typed_regions.md` §6.1); and
//! - a high-level **region seam** — `RegionScanner`, a streaming iterator that yields
//!   a gap-free repeat/satellite/unique tiling of the reference, built on the above.
//!
//! The scanner holds **no consumer policy** (no purity floor, homopolymer rule, or
//! period ceiling): a consumer passes the period range it wants, and that plus the two
//! scoring weights is the whole configuration surface.
//!
//! **Clean-room note:** derived only from Benson (1999)'s published idea (compare to a
//! period-shifted copy; grow runs of matches) and our own code — **not** from the
//! AGPL-v3 `TRF-mod` source, which was not read.
//!
//! Build status (incremental, per the impl plan): **the scanner plan's Milestones A–D are
//! done** — the type vocabulary, the `find_tandem_repeats` interval finder (lag-`p`
//! scoring + a Ruzzo–Tompa maximal-scoring-segment pass), the `RegionScanner` region seam
//! (coverage-merge tiling, resident via `over_slice` and windowed-streaming via `stream`),
//! and trf-mod parity. **Plus typed-regions B1**: `scan_windowed` promoted from a private
//! whole-contig-eager helper to a public per-window stream (`typed_regions.md` §6.1).
//!
//! Visibility note: types are `pub` (matching the sibling ng modules, e.g. `types.rs`),
//! not the arch doc's illustrative `pub(crate)` — the ng convention exposes its module
//! vocabulary as reachable API, which also lets Milestone-A types be defined ahead of
//! their Milestone-B/C consumers without tripping `dead_code`.

use crate::fasta::ChromRefFetchError;

// ---------------------------------------------------------------------
// Inputs — the scope and scoring knobs
// ---------------------------------------------------------------------

/// The inclusive range of repeat **periods** (motif lengths, in bp) to scan for.
///
/// The caller chooses the range; the algorithm imposes no ceiling of its own. This is
/// the one **constrained** type in the module — a period of 0 is meaningless and an
/// empty range (`min > max`) is a caller bug — so its fields are private and every
/// value is built through the checked [`PeriodRange::new`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PeriodRange {
    min: u8,
    max: u8,
}

impl PeriodRange {
    /// Build a period range, rejecting `min == 0` and `min > max`. These are caller
    /// bugs (period 0 has no meaning; an empty range scans nothing), so they fail
    /// loudly rather than silently coercing.
    pub fn new(min: u8, max: u8) -> Result<Self, PeriodRangeError> {
        if min == 0 {
            return Err(PeriodRangeError::ZeroMin);
        }
        if min > max {
            return Err(PeriodRangeError::MinExceedsMax { min, max });
        }
        Ok(Self { min, max })
    }

    /// The smallest period scanned (inclusive).
    #[inline]
    pub fn min(self) -> u8 {
        self.min
    }

    /// The largest period scanned (inclusive).
    #[inline]
    pub fn max(self) -> u8 {
        self.max
    }
}

/// A [`PeriodRange`] could not be built — a caller bug, never a data condition.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PeriodRangeError {
    /// `min` was 0; period 0 is meaningless.
    #[error("period range min must be >= 1 (period 0 is meaningless)")]
    ZeroMin,
    /// `min` exceeded `max`; the range would scan nothing.
    #[error("period range min ({min}) exceeds max ({max})")]
    MinExceedsMax { min: u8, max: u8 },
}

/// The scanner's general default period range, 1..=6 — a named constant, not a magic
/// literal. Each consumer overrides as it sees fit (the catalog / ng snp-str caller
/// uses 2..=6); the value may change after empirical tuning.
pub const DEFAULT_PERIODS: (u8, u8) = (1, 6);

/// Default match reward: a lag-`p` match adds this to the running score.
pub const DEFAULT_MATCH_REWARD: i32 = 2;
/// Default mismatch penalty: a lag-`p` mismatch (or a non-ACGT base) subtracts this.
/// The ratio `mismatch_penalty / match_reward` is the mismatch tolerance — a ratio `r`
/// holds a tract to purity `r / (1 + r)`; 7/2 = 3.5 holds tracts to ≈ 0.78 (spec §3.3).
/// A starting point, tuned from experiments.
pub const DEFAULT_MISMATCH_PENALTY: i32 = 7;
/// Default minimum copy count for a segment to be emitted — a repeat needs at least
/// two copies to exist at all.
pub const DEFAULT_MIN_COPIES: u32 = 2;

/// The lag-`p` self-comparison scoring plus the minimum-copies emission floor.
///
/// A tandem repeat of period `p` shows up as a high-scoring run of the signal that adds
/// [`match_reward`](Self::match_reward) where `seq[j] == seq[j-p]` (both ACGT,
/// case-insensitive) and subtracts [`mismatch_penalty`](Self::mismatch_penalty)
/// otherwise. `Default` is the scanner's general starting weights (`2 / 7 / 2`), which a
/// consumer may override with weights matched to its own purity target.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanParams {
    /// Score added on a lag-`p` match. Expected `> 0`.
    pub match_reward: i32,
    /// Score subtracted on a lag-`p` mismatch (or non-ACGT base). Expected `> 0`.
    pub mismatch_penalty: i32,
    /// A segment shorter than `min_copies · period` is not emitted (a permissive floor;
    /// a consumer applies its own, stricter, per-period floors afterward).
    ///
    /// `u32`, not `u64`: this is a **count**, not a coordinate or a length, so
    /// spec §4's widening does not reach it (the same reason ids stay `u32`, and
    /// matching `region_typing::segment_criteria::MinCopies`).
    pub min_copies: u32,
}

impl Default for ScanParams {
    fn default() -> Self {
        Self {
            match_reward: DEFAULT_MATCH_REWARD,
            mismatch_penalty: DEFAULT_MISMATCH_PENALTY,
            min_copies: DEFAULT_MIN_COPIES,
        }
    }
}

/// Above this merged-repeat length a region is [`Region::Satellite`], not a genotypeable
/// STR (satellite DNA; spec §3.6). Also bounds the streaming window halo. 1 kb is
/// comfortably above any real STR locus.
pub const DEFAULT_MAX_REPEAT_LEN: u64 = 1000;
/// Default streaming window core (bp): a memory/I/O-granularity knob only — the region
/// tiling is invariant to it.
pub const DEFAULT_WINDOW_BP: u64 = 100_000;
/// Default unique-gap bridging (0 = off): smoothing for the region tiling.
pub const DEFAULT_MERGE_GAP: u64 = 0;
/// Default minimum repeat-region length (0 = off): smoothing for the region tiling.
pub const DEFAULT_MIN_REPEAT_LEN: u64 = 0;

/// Region-seam shaping and streaming knobs — separate from [`ScanParams`] because they
/// shape the **partition and the walk**, not the detection.
///
/// `Default` = `{ 1000, 100_000, 0, 0 }`: the 1 kb satellite cap, a 100 kb streaming
/// window, and both smoothing knobs off (a pure coverage partition).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SegmentOptions {
    /// Merged repeat coverage longer than this is emitted as [`Region::Satellite`].
    pub max_repeat_len: u64,
    /// The streaming window core in bp; a memory knob, region-count-invariant.
    pub window_bp: u64,
    /// Smoothing: bridge unique gaps shorter than this into the flanking repeat.
    pub merge_gap: u64,
    /// Smoothing: reclassify a repeat region shorter than this as unique.
    pub min_repeat_len: u64,
}

impl Default for SegmentOptions {
    fn default() -> Self {
        Self {
            max_repeat_len: DEFAULT_MAX_REPEAT_LEN,
            window_bp: DEFAULT_WINDOW_BP,
            merge_gap: DEFAULT_MERGE_GAP,
            min_repeat_len: DEFAULT_MIN_REPEAT_LEN,
        }
    }
}

// ---------------------------------------------------------------------
// Output — the interval finder's natural shape
// ---------------------------------------------------------------------

/// One detected tandem-repeat interval.
///
/// Coordinates are **0-based half-open** (`[start, end)`); `period` is the lag it was
/// found at; `score` is the Ruzzo–Tompa segment total (an alignment-style
/// length-and-purity proxy). Nothing consumer-specific here — it happens to carry
/// exactly the four fields the STR post-filter reads, which is convenience, not coupling.
///
/// **0-based, deliberately, where the rest of ng is 1-based** (spec §4 says 1-based
/// *everywhere*; this is the reasoned exception, owner-agreed 2026-07-16 at B2).
/// This is a **detector/slice type**: its coordinates index the byte slice
/// `find_tandem_repeats` was handed, and `seq[iv.start..iv.end]` is what every
/// producer and consumer of it actually does. Making it 1-based would put a `- 1`
/// at every one of those sites to buy consistency with a *genetic* coordinate it
/// never expresses. The 1-based boundary belongs where the genetics starts —
/// `Locus`, `GenomeRegion` — which is exactly where §4's payoff lands. `u64` is
/// the width regardless.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepeatInterval {
    /// Tract start, inclusive (0-based).
    pub start: u64,
    /// Tract end, exclusive (0-based).
    pub end: u64,
    /// Repeat period (motif length, in bp).
    pub period: u8,
    /// Ruzzo–Tompa segment score.
    pub score: i32,
}

// ---------------------------------------------------------------------
// Output — the region seam's tiling
// ---------------------------------------------------------------------

/// A half-open, 0-based coordinate pair (`[start, end)`) — the plain span the
/// [`Region::Satellite`] / [`Region::Unique`] variants carry and the union
/// [`RepeatRegion::span`] holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RegionSpan {
    /// Start, inclusive (0-based).
    pub start: u64,
    /// End, exclusive (0-based).
    pub end: u64,
}

/// A genotypeable tandem repeat: merged repeat coverage no longer than
/// [`SegmentOptions::max_repeat_len`]. It carries the merged span **and** the
/// constituent intervals (periods may differ), so the caller can read the repeat
/// structure without re-scanning.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RepeatRegion {
    /// The union of `intervals`' spans.
    pub span: RegionSpan,
    /// The overlapping intervals merged into this region — non-empty, coordinate-ordered.
    pub intervals: Box<[RepeatInterval]>,
}

/// One tile of the reference produced by [`RegionScanner`]. The tiles a scan yields are
/// coordinate-ordered, pairwise non-overlapping, and cover the scanned span exactly;
/// consecutive tiles never share a kind (each is maximal).
///
/// The three kinds are the caller's three routes: `Repeat` → the STR path, `Unique` →
/// the SNP path, and `Satellite` → neither (mask/skip).
///
/// (`RegionScanner` lands in Milestone C; this enum is its yielded item, defined here
/// with the rest of the vocabulary.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Region {
    /// A genotypeable tandem repeat (coverage ≤ `max_repeat_len`).
    Repeat(RepeatRegion),
    /// Repeat coverage longer than `max_repeat_len` — satellite DNA. Repetitive (so not
    /// unique sequence) but too long to be a genotypeable STR, hence excluded from STR
    /// analysis.
    Satellite(RegionSpan),
    /// No tandem-repeat coverage — unique sequence.
    Unique(RegionSpan),
}

/// A region-seam iteration error. Only a reference-read failure surfaces here — the
/// detection itself is infallible; a mid-scan [`ChromRefFetchError`] (corrupt/truncated
/// reference) is fatal to the scan and fuses the iterator (spec §8). Construction-time
/// validation of the period range is [`PeriodRangeError`], separately.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum ScanError {
    /// Reading the reference failed mid-scan.
    #[error("reading the reference failed mid-scan")]
    Fetch {
        #[source]
        source: ChromRefFetchError,
    },
}

// ---------------------------------------------------------------------
// The interval finder (Milestone B) — lag-p scoring + Ruzzo–Tompa
// ---------------------------------------------------------------------

/// Canonicalise an ACGT base to an uppercase byte, or `None` for any non-ACGT byte.
/// A tandem-repeat match requires **both** compared bases to canonicalise equal (so an
/// `N == N` never matches — the two `None`s are not equal by this function's use).
#[inline]
/// Whether `motif` is **primitive** — it does not tile under any shorter period that
/// divides its length. `AAA` is not (it tiles under `A`); `ATAT` is not (tiles under
/// `AT`); `AT` and `CAG` are. A period-1 motif is always primitive. Case-insensitive,
/// matching the scan's own base comparison, so a soft-masked homopolymer still reduces.
///
/// This is what makes [`find_tandem_repeats`] honour its period range on each tract's
/// *true* period: a non-primitive motif is an alias of a shorter-period repeat, and the
/// shorter one is what the caller asked about.
fn is_primitive_motif(motif: &[u8]) -> bool {
    let p = motif.len();
    for q in 1..p {
        if !p.is_multiple_of(q) {
            continue; // only proper *divisors* of p can tile it exactly
        }
        if (q..p).all(|i| CANONICAL[motif[i] as usize] == CANONICAL[motif[i - q] as usize]) {
            return false; // tiles under its first q bases → not primitive
        }
    }
    true
}

/// Canonical base for every byte: `A`/`C`/`G`/`T` in either case map to the uppercase
/// base, and **everything else maps to `0`** — a sentinel no canonical base can equal.
///
/// This replaced a `match` (2026-07-20). The match was **54.7% of the walk's remaining
/// self-time** once the two quadratics were gone, which is what a chain of comparisons
/// costs when it runs twice per scored position: ~9.4 billion calls over tomato
/// (2 × 6 periods × 782 Mb). A 256-byte table is one L1-resident load instead.
///
/// The `0` sentinel is what lets the table serve both callers, whose semantics differ:
///
/// - [`is_primitive_motif`] compared `Option`s, so two non-ACGT bytes compared **equal**
///   (`None == None`). `0 == 0` preserves that.
/// - the scoring closure in [`find_tandem_repeats`] required `(Some(a), Some(b))` with
///   `a == b`, so a non-ACGT byte matched **nothing**, not even another one. The
///   explicit `!= 0` guard in [`bases_match`] preserves that.
static CANONICAL: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'A';
    table[b'a' as usize] = b'A';
    table[b'C' as usize] = b'C';
    table[b'c' as usize] = b'C';
    table[b'G' as usize] = b'G';
    table[b'g' as usize] = b'G';
    table[b'T' as usize] = b'T';
    table[b't' as usize] = b'T';
    table
};

/// Do two bytes canonicalise to the **same ACGT base**? A non-ACGT byte matches nothing
/// — an `N` beside an `N` is absence of evidence, not evidence of a repeat.
#[inline]
#[allow(dead_code)]
fn bases_match(x: u8, y: u8) -> bool {
    let cx = CANONICAL[x as usize];
    cx != 0 && cx == CANONICAL[y as usize]
}

/// One open segment in the Ruzzo–Tompa working stack: `l`/`r` are the cumulative score
/// strictly before the segment and through its end (so the segment's own score is
/// `r - l`); `start`/`end` are inclusive indices into the score sequence.
struct RtSeg {
    l: i64,
    r: i64,
    start: usize,
    end: usize,
    /// **Previous-smaller-element link**: `1 +` the index of the nearest segment *below*
    /// this one whose `l` is smaller, or `0` when there is none. The `+ 1` bias keeps
    /// "none" representable without an `Option`, and makes the search loop's bound check
    /// and its jump the same comparison.
    ///
    /// This is what lets the search skip rather than scan — see the search loop in
    /// [`maximal_scoring_subsequences`]. It is free to maintain: the link a segment needs
    /// is exactly the `j` the search that placed it already found.
    prev_smaller: usize,
}

/// Find every **maximal scoring subsequence** of a score sequence (Ruzzo & Tompa 1999),
/// calling `emit(start, end, score)` once per segment — `start`/`end` inclusive indices
/// into the input, `score = r - l > 0`. Segments are emitted in ascending `start` order.
///
/// A maximal scoring subsequence is a contiguous run of positive total that cannot be
/// extended or trimmed without lowering its score; substitutions and short interruptions
/// stay inside one segment as long as the surrounding score outweighs them (spec §3.4).
///
/// This is the textbook offline pass: a working stack of open segments, drained at the end.
/// Its correctness is pinned by a brute-force property test (`ruzzo_tompa_matches_brute_force`).
/// The stack holds only the segments opened since the last Rule-2 push, which is when every
/// segment below the pushed one becomes final and is emitted (see that arm).
fn maximal_scoring_subsequences(
    scores: impl Iterator<Item = i64>,
    mut emit: impl FnMut(usize, usize, i64),
) {
    let mut stack: Vec<RtSeg> = Vec::new();
    let mut total: i64 = 0;
    for (i, s) in scores.enumerate() {
        let l = total;
        total += s;
        let r = total;
        if s <= 0 {
            continue; // non-positive scores never start a segment; they advance the total
        }
        let mut cur = RtSeg {
            l,
            r,
            start: i,
            end: i,
            prev_smaller: 0,
        };
        loop {
            // The rightmost stacked segment whose `l` is smaller than `cur.l`, found by
            // **hopping down the previous-smaller chain** instead of scanning.
            //
            // Correctness: the chain skips exactly the segments that cannot be the
            // answer. Everything between `t` and `stack[t].prev_smaller` has `l >= l[t]`,
            // by that link's definition — so once `l[t] >= cur.l`, none of them is below
            // `cur.l` either, and the first `t` reached with `l[t] < cur.l` is the
            // rightmost such segment. Same answer as the scan it replaces.
            //
            // **Not a binary search** (tried, reverted, 2026-07-20): the stack's `l`
            // values are *not* sorted, and Rule 2 is what unsorts them — when nothing is
            // smaller, `cur` is pushed anyway, landing a smaller `l` above larger ones.
            // `partition_point` there compiles, passes a `debug_assert` that is compiled
            // out of release, and silently finds **zero** loci on real sequence.
            //
            // Why it was worth doing: the scan measured **139,995,860,286 steps, 7,634
            // per search, stack depth to 21,683** (tomato 8 Mb fixture, 2026-07-20) and
            // was ~100% of the walk. The costly case is the Rule-2 miss, which walked the
            // whole stack to conclude nothing qualified — and a descending stack is the
            // common shape, since mismatches drift the running total down. That case is
            // now **one hop**.
            let mut idx = stack.len();
            let mut found = None;
            while idx > 0 {
                let t = idx - 1;
                if stack[t].l < cur.l {
                    found = Some(t);
                    break;
                }
                idx = stack[t].prev_smaller;
            }
            match found {
                None => {
                    // Rule 2: no left-smaller segment → `cur` is separate. **And everything
                    // below it is final**: every stacked `l` is at least `cur.l`, so a later
                    // search that would find one of them finds `cur`, which sits above them,
                    // first — and `cur` being found, or merged into, leaves them untouched.
                    // Emitting them now keeps the stack to the segments since the last such
                    // push instead of a whole contig's: on tomato's chromosome 1 it reached
                    // 1.34 GB at period 1 when they waited for the end. They come out in the
                    // order the final drain would have emitted them.
                    for seg in &stack {
                        emit(seg.start, seg.end, seg.r - seg.l);
                    }
                    stack.clear();
                    cur.prev_smaller = 0;
                    break;
                }
                Some(j) => {
                    if stack[j].r >= cur.r {
                        // Rule 3: `cur` stays a separate segment to the right — and `j`,
                        // the answer this search just produced, is exactly the link `cur`
                        // must carry, so maintaining the chain costs nothing.
                        cur.prev_smaller = j + 1;
                        break;
                    }
                    // Rule 4: `cur` absorbs `stack[j..]` — extend its left edge and
                    // re-search. Truncating cannot invalidate a surviving link: every
                    // link points strictly below its own segment, so all of them stay
                    // within the retained prefix.
                    cur.l = stack[j].l;
                    cur.start = stack[j].start;
                    stack.truncate(j);
                }
            }
        }
        stack.push(cur);
    }
    for seg in stack.drain(..) {
        emit(seg.start, seg.end, seg.r - seg.l);
    }
}

/// Split one maximal scoring segment wherever the **motif changes**, emitting each piece as
/// `(k0, k1, score)` — the same score-index form [`maximal_scoring_subsequences`] emits.
///
/// The lag-`p` comparison asks only "is this base the one `p` back?", which is as true inside a
/// run of `A` as inside a run of `T`. So `A×12 T×12` scores 11 matches, one mismatch, 11 matches
/// at period 1, and the segment runs across the junction: one tract, motif `A`, half of it `T`.
/// The consumer's purity floor then rejects it and **both** repeats are lost. The same holds
/// for `(CA)×10 (GT)×10` at period 2, or `(GA)×12` running into `C×12`.
///
/// The test is what tells a motif change from the interruptions §3.4 keeps whole. A
/// **solid block** is a run of consecutive lag-`p` matches long enough to be a repeat by the
/// scanner's own floor — at least `min_copies` copies, and never less than two, so its motif is
/// fully determined. Its motif is compared up to rotation:
///
/// - a substitution leaves fewer than `p` matches between its two mismatches (spec §3.3), so it
///   never makes a solid block, and the blocks either side share the motif;
/// - an indel shifts phase, so the block after it carries a **rotation** of the motif before —
///   the same motif;
/// - a different repeat abutting this one carries a motif that is no rotation of it — a split.
///
/// Consecutive solid blocks with the same motif form one piece; the bases between two pieces
/// (the junction's mismatches) belong to neither.
///
/// **A short stretch of another motif inside a repeat is an interruption, not a second
/// repeat.** `ATATATATACACATATATAT` is an `AT` repeat with two copies changed to `AC`; the
/// `AC` pair is a solid block, but splitting there leaves pieces too short for a copy floor,
/// and the repeat is lost — on tomato, 876 STR segments of this shape. So when a piece's motif
/// comes back, the pieces between are **absorbed**: the nearest earlier piece of the same motif,
/// everything after it, and the new piece become one, provided the joined tract, compared base
/// by base with that motif tiled in phase, is at least as pure as the scoring holds a tract to
/// (spec §3.3: a ratio `r` holds tracts to purity `r / (1 + r)`, 7/9 ≈ 0.78 at the defaults).
/// That example matches `AT` at 18 of 20 bases, 0.9, and stays whole; so does
/// `(TAC)×2 (TGC)×4 (TAC)×3`, 0.85, whose inner motif differs from the outer by one base; so
/// does `(AAT)×4 (GAT)×3 (GGT)×2 (AAT)×4`, 0.82, with two different interruptions in a row.
/// Base by base, not piece by piece, because that is how the typing measures purity
/// afterwards. No new knob: the tolerance is the one the scoring already states.
///
/// The phase is the one at the outer piece's **last** solid block, so an indel earlier in that
/// piece makes the bases before it count as impure, and the interruption after it may then not
/// be absorbed. That is deliberate rather than overlooked: the typing measures purity in one
/// phase from the tract's start, so a joined tract with an indel in it would be rejected there
/// whole, as it was before this split existed; splitting at least keeps the pure pieces.
///
/// **A piece at either end is never absorbed**, however short: nothing of the first motif
/// follows it, so it is where one repeat stops and another starts. `G×60 A×15` is `G` over 0.8
/// of its bases, and still two repeats.
///
/// A segment with one motif throughout — every segment that is not a junction — and a segment
/// whose pieces were all absorbed are emitted unchanged, score and all. The caller sends only
/// segments with a mismatch (one without cannot change motif) that are long enough to keep.
///
/// **Out of line on purpose**: inlined into the scan's per-segment closure it made that closure
/// too big to inline into the scan loop, and the whole build paid a call per segment. Measured
/// on tomato's chromosome 1: 1.69 s and 22.5 G instructions out of line, 1.77 s and 24.6 G
/// inlined, four interleaved runs each.
#[inline(never)]
fn split_at_motif_changes(
    canonical: &[u8],
    p: usize,
    k0: usize,
    k1: usize,
    score: i64,
    params: &ScanParams,
    emit: &mut impl FnMut(usize, usize, i64),
) {
    let is_match = |k: usize| {
        let h = canonical[k + p];
        h != 0 && h == canonical[k]
    };
    // A block of `len` consecutive matches covers `len + p` bases, so `min_copies` copies need
    // `(min_copies - 1) · p` matches; at least `p`, so the block holds two whole copies.
    let solid_len = p * (params.min_copies.max(2) as usize - 1);
    let same_motif = |a: usize, b: usize| {
        let (x, y) = (&canonical[a..a + p], &canonical[b..b + p]);
        (0..p).any(|r| (0..p).all(|i| x[i] == y[(i + r) % p]))
    };

    let reward = i64::from(params.match_reward);
    let penalty = i64::from(params.mismatch_penalty);
    // How many bases in `from..to` match the motif read at `motif_at`, tiled in that block's
    // phase — before it as well as after, so a piece that opens on a fragment shorter than a
    // solid block is measured against the motif it actually repeats.
    let in_phase = |motif_at: usize, from: usize, to: usize| {
        (from..to)
            .filter(|&x| canonical[x] == canonical[motif_at + (x % p + p - motif_at % p) % p])
            .count() as i64
    };
    // Every piece as split, `(first score index, last score index)`: the second reading of a
    // joined tract. Empty — and never allocated — for a segment with one motif throughout.
    let mut split: Vec<(usize, usize)> = Vec::new();
    // The pieces so far after absorption, left to right, as `(first score index, last score
    // index, where its motif is read, bases of its tract matching that motif in phase, the first
    // of `split` it holds)`. A piece's tract is its score indices widened by `p + 1`, as the
    // segment's is.
    let mut pieces: Vec<(usize, usize, usize, i64, usize)> = Vec::new();
    // Push a piece, then absorb every interruption it closes: the nearest earlier piece of the
    // same motif, with everything between, when the joined tract matches that motif in phase
    // at the scoring's tolerance — matches / bases >= penalty / (reward + penalty), in
    // integers. Only the new bases are compared, so a long run of absorbed interruptions costs
    // one pass, not one per join.
    let push = |pieces: &mut Vec<(usize, usize, usize, i64, usize)>,
                split: &mut Vec<(usize, usize)>,
                (a, b, m): (usize, usize, usize)| {
        pieces.push((a, b, m, in_phase(m, a, b + p + 1), split.len()));
        split.push((a, b));
        loop {
            let right = pieces[pieces.len() - 1];
            let Some(j) = (0..pieces.len() - 1)
                .rev()
                .find(|&j| same_motif(pieces[j].2, right.2))
            else {
                return;
            };
            let left = pieces[j];
            let (left_end, right_end) = (left.1 + p + 1, right.1 + p + 1);
            let matches = left.3 + in_phase(left.2, left_end, right_end.max(left_end));
            if matches * (reward + penalty) < penalty * (right_end - left.0) as i64 {
                return;
            }
            pieces.truncate(j);
            pieces.push((left.0, right.1, left.2, matches, left.4));
        }
    };
    // The open piece: where it starts, the last solid block's end, and that block's start
    // (whose first `p` bases are the piece's current motif).
    let mut piece_start = k0;
    let mut motif_at: Option<usize> = None;
    let mut last_solid_end = k0;
    let mut k = k0;
    while k <= k1 {
        if !is_match(k) {
            k += 1;
            continue;
        }
        let b0 = k;
        while k <= k1 && is_match(k) {
            k += 1;
        }
        let b1 = k - 1;
        if b1 + 1 - b0 < solid_len {
            continue;
        }
        match motif_at {
            Some(m) if !same_motif(m, b0) => {
                push(&mut pieces, &mut split, (piece_start, last_solid_end, m));
                piece_start = b0;
            }
            _ => {}
        }
        motif_at = Some(b0);
        last_solid_end = b1;
    }
    let (Some(last_motif), false) = (motif_at, pieces.is_empty()) else {
        emit(k0, k1, score);
        return;
    };
    push(&mut pieces, &mut split, (piece_start, k1, last_motif));
    // Each piece's own score, over its own score indices, as the segment's was — the
    // segment's own for a piece that is the whole segment.
    let score_of = |a: usize, b: usize| {
        if (a, b) == (k0, k1) {
            score
        } else {
            (a..=b)
                .map(|k| if is_match(k) { reward } else { -penalty })
                .sum()
        }
    };
    // A joined piece is emitted **with the pieces it was joined from**: both readings, the
    // consumer's purity floor to choose between them (`segment_criteria::choose_readings`).
    for (i, &(a, b, _, _, first_split)) in pieces.iter().enumerate() {
        emit(a, b, score_of(a, b));
        let end_split = pieces.get(i + 1).map_or(split.len(), |next| next.4);
        if end_split - first_split > 1 {
            for &(sa, sb) in &split[first_split..end_split] {
                emit(sa, sb, score_of(sa, sb));
            }
        }
    }
}

/// Find every tandem-repeat interval in `seq` whose period lies in `periods` (arch §2.1).
///
/// For each period `p` it scores position `j` (`j >= p`) `+match_reward` when `seq[j]`
/// and `seq[j-p]` canonicalise to the same ACGT base and `-mismatch_penalty` otherwise
/// (a non-ACGT base never matches), takes every Ruzzo–Tompa maximal scoring segment, maps
/// it to the 0-based half-open tract it certifies, and emits it when its implied copy
/// count clears `params.min_copies`. Pure and total over arbitrary bytes; deterministic.
///
/// **Each repeat is reported once, at its true period** — the scanner honours its
/// interface, and both mechanisms it uses are pure geometry (no significance policy):
/// - a **non-primitive motif** is dropped at emission (`AA` reduces to `A`), so a
///   homopolymer is period 1, never a period-2 alias of itself;
/// - a **period-multiple re-detection of the same tract** is dropped in a post-pass:
///   a `(GAA)` period-3 repeat also matches, sloppily, as a period-6 `AAGATG` — a
///   primitive 6-mer the motif check keeps, but the period-3 detection *covers most of
///   it*, so it is the same tract seen worse.
///
/// The "covers most" test is what keeps this out of policy: a genuine shorter repeat
/// that merely *sits inside* a longer one (the `AAA` inside `(GAAA)*n`) covers only a
/// sliver of it, so it is **not** treated as a re-detection — both are kept, and which
/// is significant is the consumer's copy-floor call, not the scanner's.
pub fn find_tandem_repeats(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
) -> Vec<RepeatInterval> {
    find_tandem_repeats_for_copy_floor(seq, periods, params, |_| 0)
}

/// [`find_tandem_repeats`] for a consumer that will drop every tract holding fewer than
/// `copy_floor(period)` copies: every tract it returns at or above that floor is exactly one
/// `find_tandem_repeats` returns, and of the tracts below it, it leaves out those too short
/// to matter.
///
/// **Not simply the floor, because a tract below it can still remove one above it.** The
/// period-multiple pass drops a tract that a shorter-period tract covers half of, so a
/// period-2 tract of 9 bases, itself below a floor of 5 copies, can remove a period-4 tract
/// of 16 that clears a floor of 4. What is left out is only what cannot do that: a tract of
/// period `p` is needed when it clears its own floor, or when it is at least half as long
/// as the shortest needed tract at some multiple of `p` in the range — and those lengths are
/// worked out from the longest period down. A tract shorter than that eliminates nothing
/// that is needed, and the tracts it would have shielded are shorter still, so leaving it
/// out changes nothing that clears the floor.
///
/// Why it exists: on tomato's chromosome 1 the catalog keeps 774,585 of the 23.6 million
/// tracts `find_tandem_repeats` returns, and holding all of them through the period-multiple
/// pass was most of the build's peak memory.
pub fn find_tandem_repeats_for_copy_floor(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
    copy_floor: impl Fn(u8) -> u32,
) -> Vec<RepeatInterval> {
    // `needed[p]`: the shortest tract of period `p` that can matter, as above.
    let mut needed = vec![0u64; usize::from(periods.max()) + 1];
    for p in (periods.min()..=periods.max()).rev() {
        let own = u64::from(copy_floor(p)) * u64::from(p);
        needed[usize::from(p)] = (2 * usize::from(p)..=usize::from(periods.max()))
            .step_by(usize::from(p))
            .map(|m| needed[m].div_ceil(2))
            .fold(own, u64::min);
    }

    let mut out = Vec::new();
    let n = seq.len();
    // Canonicalise ONCE, not twice per period. `bases_match` did two dependent table
    // loads per scored position; over 6 periods that is 12 canonicalisations of every
    // base, all of the same byte.
    let mut canonical: Vec<u8> = Vec::with_capacity(n);
    canonical.extend(seq.iter().map(|&b| CANONICAL[b as usize]));
    for period in periods.min()..=periods.max() {
        let p = period as usize;
        if p >= n {
            continue; // no position has a partner `p` bases back
        }
        // Score index `k` (0-based over `p..n`) corresponds to position `j = k + p`.
        let reward = i64::from(params.match_reward);
        let penalty = -i64::from(params.mismatch_penalty);
        let back = &canonical[..n - p];
        let here = &canonical[p..];
        let scores = back
            .iter()
            .zip(here.iter())
            .map(|(&b, &h)| if h != 0 && h == b { reward } else { penalty });
        let mut emit = |k0: usize, k1: usize, score: i64| {
            // Segment [k0, k1] → tract [k0, k1 + p + 1): the earliest base involved is
            // `j0 - p = k0`, the latest is `j1 = k1 + p`, so the exclusive end is `k1+p+1`.
            let tract_start = k0 as u64;
            let tract_end = (k1 + p + 1) as u64;
            let copies = (tract_end - tract_start) / u64::from(period);
            // Emit only PRIMITIVE-period repeats — honour the requested period range on
            // the tract's true period, not on a non-primitive alias. `AAAAAA` tiles under
            // `AA`, so it is *detected* at period 2, but its motif `AA` reduces to `A`
            // (period 1): it is a homopolymer, not a dinucleotide, and reporting it as
            // period 2 would be the scanner failing its own contract. The caller asked for
            // period 2; a period-1 tract dressed as period 2 is not that.
            let motif = &seq[tract_start as usize..tract_start as usize + p];
            if copies >= u64::from(params.min_copies)
                && tract_end - tract_start >= needed[p]
                && is_primitive_motif(motif)
            {
                out.push(RepeatInterval {
                    start: tract_start,
                    end: tract_end,
                    period,
                    score: i32::try_from(score).unwrap_or(i32::MAX),
                });
            }
        };
        // The shortest tract `emit` keeps at this period. Every piece of a split lies inside its
        // segment, so a segment shorter than this yields nothing however it splits — and on
        // tomato's chromosome 1, 84 million of the 92 million segments at periods 2 to 6 are.
        let shortest = (u64::from(params.min_copies) * p as u64).max(needed[p]) as usize;
        let reward = i64::from(params.match_reward);
        maximal_scoring_subsequences(scores, |k0, k1, score| {
            if k1 + p + 1 - k0 < shortest {
                return;
            }
            // No mismatch, no motif change: a segment's score is all rewards exactly when
            // every comparison in it matched, and it is emitted whole without the walk.
            if score == reward * (k1 - k0 + 1) as i64 {
                emit(k0, k1, score);
            } else {
                split_at_motif_changes(&canonical, p, k0, k1, score, params, &mut emit);
            }
        });
    }

    // Drop a period-multiple re-detection of the same tract. `iv` (the longer period) is
    // redundant iff a shorter period `a` divides it and **covers at least half of it** —
    // "covers most" is what says *same tract*, not merely *overlapping*. A phase-shifted
    // period-6 `AAGATG` over a period-3 `GAA` is covered ~90% and goes; a period-1 `AAA`
    // sitting inside a period-4 `(GAAA)` covers a sliver and stays (a genuinely different,
    // smaller repeat — the copy floor, downstream, decides if it is significant). No
    // copy-number policy here: purely how much of `iv` the shorter `a` explains.
    out.sort_by_key(|iv| (iv.period, iv.start, std::cmp::Reverse(iv.end)));

    // Survivors grouped by period, which is what turns the redundancy test from a
    // rescan of every survivor into a binary search.
    //
    // Three facts make the grouping exact rather than a heuristic:
    //
    // 1. **Only a shorter, dividing period can eliminate `iv`**, so the search space
    //    is `divisors(period) \ {period}` — at most `{1, 2, 3}`, for period 6.
    // 2. **Those lists are already complete** when `iv` is reached: `out` is sorted by
    //    period first, so every shorter period was processed before it. (The old code
    //    relied on the same ordering, implicitly.)
    // 3. **Within one period the searched tracts are start- AND end-sorted.** Ruzzo–Tompa's
    //    segments are disjoint by construction, and each tract extends its segment by
    //    the same `p + 1`, so both endpoints stay monotonic. That is the precondition
    //    `partition_point` needs. The one exception is a joined tract's pieces, emitted
    //    beside it (`split_at_motif_changes`) and nested inside it: they are kept out of
    //    the searched lists, which loses nothing, since the tract enclosing a piece
    //    overlaps anything at least as much as the piece does. The sort puts the
    //    enclosing tract first, so a piece is known by ending no later than the reach so far.
    //
    // Why it was worth doing: the old `kept.iter().any(…)` was O(survivors²) per
    // window, and measured at **28.6% of the walk's total CPU while discarding 707 of
    // 2,084,637 intervals** (tomato 8 Mb fixture, 2026-07-20) — a near-no-op pass
    // paying a quadratic price, with up to 26,836 intervals in a single window.
    let max_period = out.iter().map(|iv| iv.period).max().unwrap_or(0) as usize;
    // Only a period with a multiple in range is ever searched (fact 1), so only those keep a
    // list, and the search reads nothing but the span. Each list is sized once to its period's
    // count in `out` (sorted by period, so two binary searches give it): grown by doubling, a
    // list of 18 M spans leaves its outgrown buffers behind as resident memory.
    let mut by_period: Vec<Vec<RegionSpan>> = (0..=max_period / 2)
        .map(|period| {
            let from = out.partition_point(|iv| usize::from(iv.period) < period);
            let to = out.partition_point(|iv| usize::from(iv.period) <= period);
            Vec::with_capacity(to - from)
        })
        .collect();
    // The furthest end of this period's searched tracts so far, and the period it belongs to.
    // Only a kept tract moves it: a joined tract dropped here leaves its pieces searchable.
    let mut reach = (0u8, 0u64);
    // Survivors are compacted in place: `retain` visits in order, once each, which is all
    // the pass needs — a second vector of survivors would hold the whole output twice.
    out.retain(|&iv| {
        let nested = reach.0 == iv.period && iv.end <= reach.1;
        let iv_len = iv.end - iv.start;
        let mut redundant = false;
        for shorter in 1..iv.period {
            if !iv.period.is_multiple_of(shorter) {
                continue;
            }
            let candidates = &by_period[shorter as usize];
            // The tracts overlapping `iv` are a contiguous run: skip those ending at
            // or before it starts, stop at the first starting at or after it ends.
            let first = candidates.partition_point(|a| a.end <= iv.start);
            for a in &candidates[first..] {
                if a.start >= iv.end {
                    break;
                }
                // `a` must cover half of `iv`, so a tract shorter than that cannot
                // qualify however it sits — cheaper than computing the overlap.
                if (a.end - a.start) * 2 < iv_len {
                    continue;
                }
                let overlap = a.end.min(iv.end).saturating_sub(a.start.max(iv.start));
                if overlap * 2 >= iv_len {
                    redundant = true;
                    break;
                }
            }
            if redundant {
                break;
            }
        }
        if !redundant && !nested {
            reach = (iv.period, iv.end);
            if let Some(list) = by_period.get_mut(usize::from(iv.period)) {
                list.push(RegionSpan {
                    start: iv.start,
                    end: iv.end,
                });
            }
        }
        !redundant
    });
    out
}

// ---------------------------------------------------------------------
// The region seam (Milestone C) — repeat / satellite / unique tiling
// ---------------------------------------------------------------------

/// A merged repeat span (the union of overlapping/abutting intervals) plus the intervals
/// that formed it — the intermediate between raw intervals and emitted [`Region`]s.
/// Union a bag of half-open `[start, end)` spans into disjoint, start-sorted spans. The
/// repeat **coverage** — the positions some repeat covers — independent of which intervals
/// produced it (so it merges identically whether fed whole-contig intervals or windowed,
/// clipped fragments).
fn union_spans(mut spans: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    spans.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (s, e) in spans {
        match out.last_mut() {
            Some(last) if s <= last.1 => last.1 = last.1.max(e),
            _ => out.push((s, e)),
        }
    }
    out
}

/// Build the final region tiling of `[0, n)` from two inputs kept **deliberately separate**:
/// the repeat `coverage` (span-level, which classifies Repeat/Satellite and fixes the
/// boundaries) and the exact `str_intervals` (which decorate each `Repeat` with its
/// constituents). Splitting them is what makes the windowed path invariant: coverage
/// coalesces satellite fragments across windows, while STR intervals — always ≤
/// `max_repeat_len`, so captured whole by one window — are attached by start.
///
/// Steps: union the coverage into spans, apply the [`SegmentOptions`] smoothing (`merge_gap`
/// bridges short unique gaps, `min_repeat_len` drops short repeat blips), classify each
/// surviving span `Repeat` vs `Satellite` by `max_repeat_len` (attaching the `str_intervals`
/// whose start lies in a `Repeat` span), and fill every gap with `Unique`. The output tiles
/// `[0, n)` exactly, coordinate-ordered, with no two same-kind tiles adjacent.
fn build_regions(
    coverage: Vec<(u64, u64)>,
    mut str_intervals: Vec<RepeatInterval>,
    n: u64,
    opts: &SegmentOptions,
) -> Vec<Region> {
    let mut spans = union_spans(coverage);

    // Smoothing — bridge unique gaps shorter than `merge_gap`.
    if opts.merge_gap > 0 {
        let mut bridged: Vec<(u64, u64)> = Vec::new();
        for (s, e) in spans {
            match bridged.last_mut() {
                Some(last) if s - last.1 < opts.merge_gap => last.1 = e,
                _ => bridged.push((s, e)),
            }
        }
        spans = bridged;
    }
    // Smoothing — drop repeat spans shorter than `min_repeat_len` (they become unique).
    if opts.min_repeat_len > 0 {
        spans.retain(|&(s, e)| e - s >= opts.min_repeat_len);
    }

    str_intervals.sort_by_key(|iv| (iv.start, iv.end));
    let mut iv_idx = 0usize; // a single forward cursor into the start-sorted intervals

    let mut out = Vec::new();
    let mut pos = 0u64;
    for (start, end) in spans {
        if start > pos {
            out.push(Region::Unique(RegionSpan {
                start: pos,
                end: start,
            }));
        }
        let region_span = RegionSpan { start, end };
        if end - start > opts.max_repeat_len {
            out.push(Region::Satellite(region_span));
        } else {
            // Attach the intervals whose start lies in this span. Skipped intervals (in
            // gaps, satellites, or dropped spans) are stepped over by the shared cursor.
            while iv_idx < str_intervals.len() && str_intervals[iv_idx].start < start {
                iv_idx += 1;
            }
            let lo = iv_idx;
            while iv_idx < str_intervals.len() && str_intervals[iv_idx].start < end {
                iv_idx += 1;
            }
            out.push(Region::Repeat(RepeatRegion {
                span: region_span,
                intervals: str_intervals[lo..iv_idx].to_vec().into_boxed_slice(),
            }));
        }
        pos = end;
    }
    if pos < n {
        out.push(Region::Unique(RegionSpan { start: pos, end: n }));
    }
    out
}

/// Tile a resident sequence into an ordered, gap-free run of [`Region`]s (spec §3.6). The
/// **oracle** the windowed streaming path (`stream`) must reproduce (window-count invariance):
/// its coverage is every interval's span and its STR intervals are all of them.
fn tile(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
    opts: &SegmentOptions,
) -> Vec<Region> {
    let n = seq.len() as u64;
    if n == 0 {
        return Vec::new();
    }
    let intervals = find_tandem_repeats(seq, periods, params);
    let coverage = intervals.iter().map(|iv| (iv.start, iv.end)).collect();
    build_regions(coverage, intervals, n, opts)
}

/// Streams a repeat/satellite/unique tiling of a reference as an iterator of [`Region`]s —
/// the routing seam the snp/ssr caller consumes (spec §3.6). Yields `Result<Region,
/// ScanError>`.
///
/// Two constructors: [`over_slice`](Self::over_slice) tiles a resident sequence, and
/// `stream` tiles a whole contig windowed through a `ChromRefFetcher`
/// without holding it in memory. `stream` is **window-count invariant for the STR tiling** —
/// it yields exactly the `Repeat`/`Unique` structure `over_slice` would, regardless of
/// `window_bp`. A `Satellite` (a repeat longer than `max_repeat_len`, hence longer than a
/// window's reach) may segment differently between the two, because its internal structure
/// depends on global context a bounded window cannot see; that is harmless, since satellites
/// are mask/skip regions, not genotyping targets (see `scan_windowed`).
pub struct RegionScanner {
    regions: std::vec::IntoIter<Region>,
}

impl RegionScanner {
    /// Tile a small, already-resident sequence (spec §3.6). Infallible (yields only `Ok`).
    pub fn over_slice(
        seq: &[u8],
        periods: PeriodRange,
        params: &ScanParams,
        opts: &SegmentOptions,
    ) -> Self {
        Self {
            regions: tile(seq, periods, params, opts).into_iter(),
        }
    }
}

impl Iterator for RegionScanner {
    type Item = Result<Region, ScanError>;

    fn next(&mut self) -> Option<Self::Item> {
        self.regions.next().map(Ok)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- PeriodRange -------------------------------------------------------

    #[test]
    fn period_range_accepts_valid_ranges() {
        assert_eq!(
            PeriodRange::new(1, 6).map(|r| (r.min(), r.max())),
            Ok((1, 6))
        );
        // A single-period range is valid.
        assert_eq!(
            PeriodRange::new(6, 6).map(|r| (r.min(), r.max())),
            Ok((6, 6))
        );
        assert_eq!(
            PeriodRange::new(2, 6).map(|r| (r.min(), r.max())),
            Ok((2, 6))
        );
    }

    #[test]
    fn period_range_rejects_zero_min() {
        assert_eq!(PeriodRange::new(0, 6), Err(PeriodRangeError::ZeroMin));
        // Zero-min is checked before the min>max rule, so (0, 0) is ZeroMin, not empty.
        assert_eq!(PeriodRange::new(0, 0), Err(PeriodRangeError::ZeroMin));
    }

    #[test]
    fn period_range_rejects_min_above_max() {
        assert_eq!(
            PeriodRange::new(3, 2),
            Err(PeriodRangeError::MinExceedsMax { min: 3, max: 2 })
        );
    }

    // ---- Defaults ----------------------------------------------------------

    #[test]
    fn default_periods_is_one_to_six() {
        assert_eq!(DEFAULT_PERIODS, (1, 6));
        // The named default range is itself a valid PeriodRange.
        assert!(PeriodRange::new(DEFAULT_PERIODS.0, DEFAULT_PERIODS.1).is_ok());
    }

    #[test]
    fn scan_params_default_is_two_seven_two() {
        let p = ScanParams::default();
        assert_eq!(p.match_reward, 2);
        assert_eq!(p.mismatch_penalty, 7);
        assert_eq!(p.min_copies, 2);
    }

    #[test]
    fn segment_options_default_is_cap_window_no_smoothing() {
        let o = SegmentOptions::default();
        assert_eq!(o.max_repeat_len, 1000);
        assert_eq!(o.window_bp, 100_000);
        assert_eq!(o.merge_gap, 0);
        assert_eq!(o.min_repeat_len, 0);
    }

    // ---- Ruzzo–Tompa vs a brute-force oracle -------------------------------

    /// Collect `maximal_scoring_subsequences` into a start-sorted `(start, end, score)`
    /// vector for comparison.
    fn rt_segments(scores: &[i64]) -> Vec<(usize, usize, i64)> {
        let mut out = Vec::new();
        maximal_scoring_subsequences(scores.iter().copied(), |a, b, s| out.push((a, b, s)));
        out.sort_unstable();
        out
    }

    /// Brute-force all maximal scoring subsequences straight from the Ruzzo–Tompa
    /// definition (O(n⁴), for small test inputs only — the independent oracle). An
    /// inclusive range `[i, j]` qualifies iff: (1) its score is positive; (2) every
    /// *proper* sub-range scores strictly less (so it is internally optimal); and (3) it
    /// is not properly contained in another range satisfying (1)+(2) (so it is maximal).
    fn brute_segments(scores: &[i64]) -> Vec<(usize, usize, i64)> {
        let n = scores.len();
        let mut pre = vec![0i64; n + 1];
        for i in 0..n {
            pre[i + 1] = pre[i] + scores[i];
        }
        let score = |i: usize, j: usize| pre[j + 1] - pre[i]; // inclusive [i, j]

        // Candidates: positive score, every proper sub-range strictly smaller.
        let mut cands: Vec<(usize, usize, i64)> = Vec::new();
        for i in 0..n {
            for j in i..n {
                let sc = score(i, j);
                if sc <= 0 {
                    continue;
                }
                let internally_optimal =
                    (i..=j).all(|i2| (i2..=j).all(|j2| (i2, j2) == (i, j) || score(i2, j2) < sc));
                if internally_optimal {
                    cands.push((i, j, sc));
                }
            }
        }

        // Keep only those not properly contained in another candidate (maximality).
        let mut res: Vec<(usize, usize, i64)> = cands
            .iter()
            .copied()
            .filter(|&(i, j, _)| {
                !cands
                    .iter()
                    .any(|&(a, b, _)| a <= i && j <= b && (a, b) != (i, j))
            })
            .collect();
        res.sort_unstable();
        res
    }

    #[test]
    fn ruzzo_tompa_matches_brute_force_on_crafted_cases() {
        let cases: &[&[i64]] = &[
            &[],
            &[-7, -7, -7],
            &[2, 2, 2, 2],
            &[2, -7, 2],  // two singletons — the interruption is too costly to bridge
            &[2, -3, 2],  // one merged segment — the dip is bridged
            &[5, -10, 8], // dip below the left L → two segments
            &[5, -3, 8],  // dip stays above → one segment
            &[2, 2, -7, 2, 2],
            &[-7, 2, 2, -7, 2, -7, 2, 2, -7],
            &[1, 1, 1, -2, 1, 1, 1, -10, 5, 5],
        ];
        for case in cases {
            assert_eq!(rt_segments(case), brute_segments(case), "case {case:?}");
        }
    }

    #[test]
    fn ruzzo_tompa_matches_brute_force_on_pseudo_random() {
        // Deterministic LCG over ±-style scores; no rand dependency.
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        for _ in 0..2000 {
            let len = (next() % 160) as usize;
            let scores: Vec<i64> = (0..len)
                .map(|_| if next() % 4 == 0 { 2 } else { -7 }) // 25% match-like
                .collect();
            assert_eq!(
                rt_segments(&scores),
                brute_segments(&scores),
                "scores {scores:?}"
            );
        }
    }

    /// What a consumer with a copy floor keeps is the same from the floor-aware scan as from
    /// the whole one — on repeat-dense sequence, where period-multiple eliminations are
    /// common, and over floors from permissive to strict.
    #[test]
    fn the_floor_aware_scan_keeps_what_the_floor_keeps() {
        let mut state: u64 = 0x2545F4914F6CDD1D;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        let floors: [[u32; 6]; 5] = [
            [5, 5, 4, 4, 4, 3],
            [2, 2, 2, 2, 2, 2],
            [3, 7, 2, 5, 2, 2],
            [12, 2, 6, 2, 3, 2],
            [9, 9, 9, 9, 9, 9],
        ];
        let mut sub_floor_kept = 0usize;
        for _ in 0..300 {
            let mut seq = Vec::new();
            while seq.len() < 3_000 {
                let period = 1 + (next() % 6) as usize;
                let motif: Vec<u8> = (0..period)
                    .map(|_| b"ACGT"[(next() % 4) as usize])
                    .collect();
                for _ in 0..1 + next() % 10 {
                    for &b in &motif {
                        seq.push(if next() % 12 == 0 {
                            b"ACGT"[(next() % 4) as usize]
                        } else {
                            b
                        });
                    }
                }
            }
            let whole = find_tandem_repeats(&seq, p(1, 6), &ScanParams::default());
            for floor in &floors {
                let at = |period: u8| floor[usize::from(period) - 1];
                let clears = |iv: &&RepeatInterval| {
                    (iv.end - iv.start) / u64::from(iv.period) >= u64::from(at(iv.period))
                };
                let want: Vec<_> = whole.iter().filter(clears).collect();
                let got_all =
                    find_tandem_repeats_for_copy_floor(&seq, p(1, 6), &ScanParams::default(), at);
                let got: Vec<_> = got_all.iter().filter(clears).collect();
                assert_eq!(got, want, "floor {floor:?}");
                sub_floor_kept += got_all.len().saturating_sub(got.len());
            }
        }
        assert!(
            sub_floor_kept > 0,
            "the sequences exercise sub-floor tracts"
        );
    }

    // ---- find_tandem_repeats ----------------------------------------------

    fn p(min: u8, max: u8) -> PeriodRange {
        PeriodRange::new(min, max).unwrap()
    }

    /// A clean (CAG)*8 tract detected at period 3 as one exact interval.
    #[test]
    fn finds_a_clean_perfect_tract() {
        let seq = b"CAGCAGCAGCAGCAGCAGCAGCAG"; // 8 copies, 24 bp
        let got = find_tandem_repeats(seq, p(3, 3), &ScanParams::default());
        assert_eq!(
            got,
            vec![RepeatInterval {
                start: 0,
                end: 24,
                period: 3,
                score: 2 * (24 - 3), // match_reward × matches
            }]
        );
    }

    /// A single substitution inside a tract keeps it **one** interval (spec §3.4).
    #[test]
    fn a_substitution_keeps_one_interval() {
        // (CAG)*5, then one base flipped, then (CAG)*5.
        let mut seq = Vec::new();
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"CAT"); // the interrupting copy (G→T)
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        let period3: Vec<_> = got.iter().filter(|r| r.period == 3).collect();
        assert_eq!(
            period3.len(),
            1,
            "one span across the substitution: {got:?}"
        );
        assert_eq!(period3[0].start, 0);
        assert_eq!(period3[0].end, seq.len() as u64);
    }

    /// A single indel (a bounded mismatch burst) also keeps the tract **one** interval —
    /// the lag-p comparison is local, so downstream copies resume in phase (spec §3.4).
    #[test]
    fn a_single_indel_keeps_one_interval() {
        // (CAG)*6, one inserted base, (CAG)*6 — phase shifts but the burst is bounded.
        let mut seq = Vec::new();
        for _ in 0..6 {
            seq.extend_from_slice(b"CAG");
        }
        seq.push(b'T'); // a 1 bp insertion
        for _ in 0..6 {
            seq.extend_from_slice(b"CAG");
        }
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        let period3: Vec<_> = got.iter().filter(|r| r.period == 3).collect();
        assert_eq!(
            period3.len(),
            1,
            "the single indel is a bounded burst, not a split: {got:?}"
        );
        assert_eq!(period3[0].start, 0);
        assert_eq!(period3[0].end, seq.len() as u64);
    }

    /// Two abutting repeats of different motifs come back as **two** intervals, each with its
    /// own motif — not one interval in the first motif, half of it the second, which the
    /// purity floor then rejects whole. The four two-repeat shapes a simulator found merged.
    #[test]
    fn abutting_repeats_of_different_motifs_are_split() {
        /// A sequence, and the `(start, end, period)` intervals it must yield.
        type Case = (Vec<u8>, Vec<(u64, u64, u8)>);
        let rep = |unit: &[u8], n: usize| unit.repeat(n);
        let cases: Vec<Case> = vec![
            (
                [rep(b"A", 12), rep(b"T", 12)].concat(),
                vec![(0, 12, 1), (12, 24, 1)],
            ),
            (
                [rep(b"G", 60), rep(b"A", 15)].concat(),
                vec![(0, 60, 1), (60, 75, 1)],
            ),
            (
                [rep(b"GA", 12), rep(b"C", 12)].concat(),
                vec![(0, 24, 2), (24, 36, 1)],
            ),
            (
                [rep(b"CA", 10), rep(b"GT", 10)].concat(),
                vec![(0, 20, 2), (20, 40, 2)],
            ),
        ];
        for (seq, want) in cases {
            let got = find_tandem_repeats(&seq, p(1, 6), &ScanParams::default());
            let mut spans: Vec<(u64, u64, u8)> =
                got.iter().map(|r| (r.start, r.end, r.period)).collect();
            spans.sort_unstable();
            for w in &want {
                assert!(
                    spans.contains(w),
                    "{} lacks {w:?}: {spans:?}",
                    String::from_utf8_lossy(&seq)
                );
            }
            // Nothing spans the junction in the first motif's period.
            let junction = want[0].1;
            assert!(
                !got.iter().any(|r| r.start < junction - 1
                    && r.end > junction + 1
                    && r.period == want[0].2),
                "{} still has a merged interval: {spans:?}",
                String::from_utf8_lossy(&seq)
            );
        }
    }

    /// The intervals of `got`, as `(start, end)`, after asserting the first is `joined` and
    /// every other one lies inside it — the shape of a joined tract emitted with its pieces.
    fn joined_with_pieces(got: &[RepeatInterval], joined: (u64, u64)) -> Vec<(u64, u64)> {
        let spans: Vec<(u64, u64)> = got.iter().map(|r| (r.start, r.end)).collect();
        assert_eq!(spans.first(), Some(&joined), "{got:?}");
        assert!(spans.len() > 1, "the pieces come with it: {got:?}");
        assert!(
            spans[1..]
                .iter()
                .all(|&(s, e)| s >= joined.0 && e <= joined.1),
            "every piece lies inside the joined tract: {got:?}"
        );
        spans
    }

    /// Two copies of another motif **inside** a repeat are an interruption, not a second
    /// repeat: `(AT)×4 A (CA)×2 (TA)×4 T` comes back as one `AT` interval over all of it, as it
    /// did before the split existed — and with the pieces it was joined from, for the typing to
    /// fall back on if the joined tract is not pure enough for it (`choose_readings`). It
    /// matches `AT` at 18 of 20 bases, above the 7/9 the default scoring holds a tract to.
    /// Found on tomato, where splitting this shape lost 876 STR segments.
    #[test]
    fn a_short_run_of_another_motif_inside_a_repeat_is_not_split() {
        let seq = b"ATATATATACACATATATAT";
        let got = find_tandem_repeats(seq, p(2, 2), &ScanParams::default());
        assert_eq!(
            joined_with_pieces(&got, (0, 20)),
            vec![(0, 20), (0, 9), (8, 13), (12, 20)]
        );
    }

    /// Two **different** interruptions in a row are absorbed too, once the outer motif comes
    /// back: `(AAT)×4 (GAT)×3 (GGT)×2 (AAT)×4` matches `AAT` at 32 of 39 bases, 0.82. A test of
    /// neighbours alone never saw it — no piece there has the same motif on both sides — and
    /// on tomato 13 STR segments of this shape became plain sequence.
    #[test]
    fn two_different_interruptions_in_a_row_are_absorbed() {
        let seq = [
            b"AAT".repeat(4),
            b"GAT".repeat(3),
            b"GGT".repeat(2),
            b"AAT".repeat(4),
        ]
        .concat();
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        joined_with_pieces(&got, (0, 39));
    }

    /// The same interruption, too long for the tolerance: `(CA)×10 (GT)×10 (CA)×10` has `CA`
    /// over two thirds of its bases, under 7/9, so all three repeats come back.
    #[test]
    fn a_long_run_of_another_motif_inside_a_repeat_is_split() {
        let seq = [b"CA".repeat(10), b"GT".repeat(10), b"CA".repeat(10)].concat();
        let got = find_tandem_repeats(&seq, p(2, 2), &ScanParams::default());
        assert_eq!(
            got.iter().map(|r| (r.start, r.end)).collect::<Vec<_>>(),
            vec![(0, 20), (20, 40), (40, 60)],
            "{got:?}"
        );
    }

    /// A homopolymer on its own is untouched by the split.
    #[test]
    fn a_lone_homopolymer_is_one_interval() {
        let seq = b"AAAAAAAAAAAA";
        let got = find_tandem_repeats(seq, p(1, 6), &ScanParams::default());
        assert_eq!(
            got.iter()
                .map(|r| (r.start, r.end, r.period))
                .collect::<Vec<_>>(),
            vec![(0, 12, 1)]
        );
    }

    /// A run of `N` yields nothing — `N == N` is not a match (spec §3.5).
    #[test]
    fn an_n_run_emits_nothing() {
        let seq = b"NNNNNNNNNNNNNNNNNNNN";
        assert!(find_tandem_repeats(seq, p(1, 6), &ScanParams::default()).is_empty());
    }

    /// A soft-masked (lowercase) tract is still found — comparison is case-insensitive.
    #[test]
    fn a_soft_masked_tract_is_found() {
        let seq = b"cagcagcagcagcagcagcagcag"; // lowercase (CAG)*8
        let got = find_tandem_repeats(seq, p(3, 3), &ScanParams::default());
        assert_eq!(got.len(), 1);
        assert_eq!((got[0].start, got[0].end, got[0].period), (0, 24, 3));
    }

    /// An (AT)*n tract matches at period 2 **and** its multiples 4 and 6 — all emitted;
    /// the finder does not de-duplicate periods (spec §3.5).
    #[test]
    fn only_the_primitive_period_is_emitted() {
        // (AT)*10 tiles under `AT` (period 2), `ATAT` (4) and `ATATAT` (6), so the scorer
        // matches at all three. But its *primitive* period is 2, and the scanner honours
        // its contract: it emits the tract once, at period 2 — not as a period-4 or
        // period-6 alias, which are the same repeat wearing a longer motif.
        let seq = b"ATATATATATATATATATAT"; // (AT)*10, 20 bp
        let got = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        let periods: std::collections::BTreeSet<u8> = got.iter().map(|r| r.period).collect();
        assert_eq!(
            periods,
            std::collections::BTreeSet::from([2]),
            "only the primitive period 2, no 4/6 aliases (nor 3/5 non-divisors): {got:?}"
        );
    }

    /// A genuine period-4 primitive (`ACGT`) is emitted at 4, not mistaken for an alias.
    #[test]
    fn a_genuine_higher_primitive_period_is_kept() {
        let seq = b"ACGTACGTACGTACGT"; // (ACGT)*4 — primitive period 4
        let got = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        assert!(
            got.iter().any(|r| r.period == 4),
            "ACGT is a primitive period-4 motif and must survive: {got:?}"
        );
        assert!(
            got.iter().all(|r| r.period != 2 && r.period != 3),
            "it does not tile under any shorter period, so no 2/3 detection: {got:?}"
        );
    }

    /// A homopolymer is a period-1 tract and nothing else — never a di/tri/… alias.
    #[test]
    fn a_homopolymer_is_only_period_one() {
        let seq = b"AAAAAAAAAAAAAAAAAAAA"; // poly-A, 20 bp
        let got = find_tandem_repeats(seq, p(1, 6), &ScanParams::default());
        let periods: std::collections::BTreeSet<u8> = got.iter().map(|r| r.period).collect();
        assert_eq!(
            periods,
            std::collections::BTreeSet::from([1]),
            "a poly-A is period 1 only — the whole aliasing problem, gone at the source: {got:?}"
        );
        // And asked for period 2+ only, the scanner returns nothing for it — the contract.
        let got2 = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        assert!(
            got2.is_empty(),
            "asked for period 2..6, a homopolymer is not any of them: {got2:?}"
        );
    }

    /// **A phase-shifted period-multiple re-detection is dropped; a repeat that merely
    /// sits inside another is kept.** This is the geometry test — the scanner resolves
    /// same-tract redundancy by *how much* the shorter covers the longer, no copy floor.
    #[test]
    fn a_period_multiple_re_detection_is_dropped_by_geometry() {
        // A (GAA)*n repeat with a couple of stray bases before it: the scanner finds
        // the clean period-3, and also a sloppy period-6 that starts earlier (motif
        // `AAGATG`). The period-3 covers most of the period-6, so the period-6 goes.
        let seq = b"ATGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAG";
        let got = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        assert!(
            got.iter().any(|r| r.period == 3),
            "the primitive period-3 GAA survives: {got:?}"
        );
        assert!(
            !got.iter().any(|r| r.period == 6),
            "the period-6 re-detection of the same tract is dropped: {got:?}"
        );

        // A homopolymer INSIDE a higher-period motif is a different, smaller repeat that
        // covers only a sliver — it must NOT eliminate its parent. (GAAA)*n has an `AAA`
        // run in every unit; the period-4 tetranucleotide must survive them.
        let seq = b"GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"; // (GAAA)*8
        let got = find_tandem_repeats(seq, p(1, 6), &ScanParams::default());
        assert!(
            got.iter().any(|r| r.period == 4),
            "the period-4 (GAAA) survives its own internal AAA homopolymer runs: {got:?}"
        );
    }

    /// The `min_copies` floor drops sub-threshold segments (e.g. a lone coincidental match).
    #[test]
    fn the_min_copies_floor_drops_short_segments() {
        // One coincidental period-2 match embedded in non-repetitive sequence.
        let seq = b"ACGTACAGTCTGAC";
        let got = find_tandem_repeats(seq, p(2, 2), &ScanParams::default());
        assert!(
            got.iter().all(|r| (r.end - r.start) / 2 >= 2),
            "every emitted interval has >= min_copies copies: {got:?}"
        );
    }

    // ---- RegionScanner::over_slice (the region seam, C1) -------------------

    /// A region's `(start, end, kind)` for the tiling-contract check — kind is 0 Repeat,
    /// 1 Satellite, 2 Unique.
    fn triple(r: &Region) -> (u64, u64, u8) {
        match r {
            Region::Repeat(rr) => (rr.span.start, rr.span.end, 0),
            Region::Satellite(s) => (s.start, s.end, 1),
            Region::Unique(s) => (s.start, s.end, 2),
        }
    }

    fn regions_of(seq: &[u8], periods: PeriodRange, opts: &SegmentOptions) -> Vec<Region> {
        RegionScanner::over_slice(seq, periods, &ScanParams::default(), opts)
            .map(|r| r.unwrap())
            .collect()
    }

    /// Every tiling is coordinate-ordered, gap-free over `[0, n)`, non-overlapping, and has
    /// no two same-kind tiles adjacent.
    fn assert_tiles(regions: &[Region], n: u64) {
        let mut pos = 0u64;
        let mut prev_kind: Option<u8> = None;
        for r in regions {
            let (s, e, kind) = triple(r);
            assert_eq!(
                s, pos,
                "tile starts at {s}, expected {pos} (gap or overlap)"
            );
            assert!(e > s, "empty tile [{s}, {e})");
            assert_ne!(Some(kind), prev_kind, "two same-kind tiles adjacent at {s}");
            prev_kind = Some(kind);
            pos = e;
        }
        assert_eq!(pos, n, "tiling must cover exactly [0, {n})");
    }

    #[test]
    fn empty_input_yields_no_regions() {
        assert!(regions_of(b"", p(1, 6), &SegmentOptions::default()).is_empty());
    }

    #[test]
    fn a_lone_repeat_tiles_as_unique_repeat_unique() {
        // Non-repetitive flanks around a clean (CAG)*8 tract.
        let mut seq = b"TTGGA".to_vec();
        for _ in 0..8 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"GGTTA");
        let regions = regions_of(&seq, p(3, 3), &SegmentOptions::default());
        assert_tiles(&regions, seq.len() as u64);
        let kinds: Vec<u8> = regions.iter().map(|r| triple(r).2).collect();
        assert_eq!(kinds, vec![2, 0, 2], "Unique, Repeat, Unique: {regions:?}");
        // The single Repeat covers the tract.
        if let Region::Repeat(rr) = &regions[1] {
            assert_eq!((rr.span.start, rr.span.end), (5, 5 + 24));
        } else {
            panic!("middle tile is the repeat");
        }
    }

    #[test]
    fn adjacent_repeats_of_different_periods_merge_into_one_region() {
        // An (AT)*8 tract abutting a (CAG)*5 tract: two primitive repeats (period 2 and 3)
        // whose coverage abuts, so the region tiling merges them into one `Repeat` region
        // carrying both intervals. (Before the primitive-period fix this test used (AT)*10's
        // 2/4/6 aliases to exercise the merge; a real repeat is now needed, because aliases
        // no longer exist.)
        let seq = b"ATATATATATATATATCAGCAGCAGCAGCAG"; // (AT)*8 then (CAG)*5
        let regions = regions_of(seq, p(2, 6), &SegmentOptions::default());
        assert_tiles(&regions, seq.len() as u64);
        let repeats: Vec<_> = regions
            .iter()
            .filter_map(|r| match r {
                Region::Repeat(rr) => Some(rr),
                _ => None,
            })
            .collect();
        assert_eq!(repeats.len(), 1, "one merged repeat region: {regions:?}");
        let periods: std::collections::BTreeSet<u8> =
            repeats[0].intervals.iter().map(|iv| iv.period).collect();
        assert_eq!(
            periods,
            std::collections::BTreeSet::from([2, 3]),
            "the merged region carries both PRIMITIVE periods, no aliases: {:?}",
            repeats[0].intervals
        );
    }

    #[test]
    fn a_repeat_over_the_cap_is_a_satellite() {
        // The same (CAG)*8 = 24 bp tract, but with a 10 bp satellite cap → Satellite.
        let mut seq = b"TTGGA".to_vec();
        for _ in 0..8 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"GGTTA");
        let opts = SegmentOptions {
            max_repeat_len: 10,
            ..SegmentOptions::default()
        };
        let regions = regions_of(&seq, p(3, 3), &opts);
        assert_tiles(&regions, seq.len() as u64);
        assert_eq!(
            regions.iter().map(|r| triple(r).2).collect::<Vec<_>>(),
            vec![2, 1, 2],
            "Unique, Satellite, Unique: {regions:?}"
        );
    }

    #[test]
    fn merge_gap_bridges_two_nearby_repeats() {
        // Two (CAG)*5 tracts separated by a 12 bp non-repetitive gap → two intervals.
        let mut seq = Vec::new();
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"TGACGTTGACGT"); // 12 bp gap the finder won't bridge
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        // Off: two separate Repeats, a Unique gap between them.
        let off = regions_of(&seq, p(3, 3), &SegmentOptions::default());
        let repeats_off = off.iter().filter(|r| triple(r).2 == 0).count();
        assert_eq!(repeats_off, 2, "without merge_gap, two repeats: {off:?}");
        assert_tiles(&off, seq.len() as u64);
        // On: a merge_gap wider than the 12 bp gap coalesces them into one Repeat.
        let opts = SegmentOptions {
            merge_gap: 20,
            ..SegmentOptions::default()
        };
        let on = regions_of(&seq, p(3, 3), &opts);
        assert_eq!(
            on.iter().filter(|r| triple(r).2 == 0).count(),
            1,
            "with merge_gap, one bridged repeat: {on:?}"
        );
        assert_tiles(&on, seq.len() as u64);
    }

    #[test]
    fn min_repeat_len_reclassifies_a_short_blip_as_unique() {
        // An isolated (CAG)*3 = 9 bp tract, flanked.
        let mut seq = b"TTGGA".to_vec();
        for _ in 0..3 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"GGTTA");
        // Off: the 9 bp tract is a Repeat.
        let off = regions_of(&seq, p(3, 3), &SegmentOptions::default());
        assert!(
            off.iter().any(|r| triple(r).2 == 0),
            "a repeat is present: {off:?}"
        );
        // On: min_repeat_len = 20 drops it → the whole contig is one Unique tile.
        let opts = SegmentOptions {
            min_repeat_len: 20,
            ..SegmentOptions::default()
        };
        let on = regions_of(&seq, p(3, 3), &opts);
        assert_eq!(
            on.iter().map(|r| triple(r).2).collect::<Vec<_>>(),
            vec![2],
            "the short repeat is reclassified unique, gaps merged: {on:?}"
        );
        assert_tiles(&on, seq.len() as u64);
    }

    // ---- scan_windowed (the streamed primitive, B1) -------------------------

    // ---- RegionScanner::stream (windowed, C2) ------------------------------
}
