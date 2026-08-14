//! The census — what each sample writes down at a kept locus.
//!
//! **The same questions put to every sample.** A census is a set of questions asked
//! identically of a whole population, which is the one property this whole route rests on:
//! the loci are the same in every sample, so the evidence gathered at them can be compared
//! sample against sample rather than only summarised within one.
//!
//! **The fit reads nothing but this evidence**, so what is missing from it cannot be
//! recovered later: the walk visits every locus once, and a field not written then would
//! need a second traversal of the reads, which is the one thing this whole step exists to
//! avoid.
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_joint_records.md`. Types:
//! `doc/devel/ng/arch/parameter_prepass_joint_records.md`.
//!
//! **Two kinds of evidence, because the two paths observe different things.** At an ordinary
//! position the observation is *which base a read showed*; at a repeat tract it is *how long
//! a tract a read showed*. Five per-base buckets cannot express a length and a length
//! distribution cannot express which base was substituted, so the two share the selection
//! rule and nothing about their contents.
//!
//! **Three states must survive a write and a read** — a locus never walked (a bug), a locus
//! walked with no coverage (data), and a locus whose reads all matched (data) — and at a
//! repeat tract there is a fourth: reads reached the locus but none crossed the whole tract.
//! [`DepthCode`] is what makes the first expressible; [`SsrEvidence::covering_not_crossing`]
//! the fourth.
//!
//! # What is not here yet
//!
//! **The byte-level census file.** Where the evidence lives between the walk and the fit is a
//! non-goal of the spec (§1.2) — what it requires is a property, that the fit reaches every
//! sample's evidence without walking the reads again. The encoding that carries the *content*
//! is here and is tested by packing and unpacking it; the framing that would put it in a file
//! is the next unit, and it changes none of these types.

use std::collections::BTreeMap;

use md5::{Digest, Md5};

use crate::ng::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use crate::ng::parameter_estimation::generic::depth_bins::{DepthBin, DepthBinEdges};
use crate::ng::parameter_estimation::joint::loci::{
    CensusLoci, CensusLociDigest, CensusLociDigester, SelectionTerms,
};
use crate::ng::repeat_catalog::StratumCounts;
use crate::ng::types::{Bp, ContigId, GenomePosition, GenomeRegion, Position, ReadGroupId};

// ---------------------------------------------------------------------
// The generic half of the census
// ---------------------------------------------------------------------

/// What a read showed at a position: one of the four bases, or anything else.
///
/// **The four are kept apart rather than collapsed to "reference and other"**, and the
/// reason is structural: the quantity fitted at a locus is *the frequency of an allele*, so
/// the allele has to exist as a thing that can have a frequency, and "non-reference" is not
/// one. It is also what lets the fit **sum over** which non-reference base is segregating
/// instead of picking the largest — picking is conditioning on a maximum, small per site,
/// one-directional, and landing on exactly the rare-frequency classes everything downstream
/// reads.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum ObservedAllele {
    A,
    C,
    G,
    T,
    /// An indel, a spanning deletion, an `N` — anything that is not one of the four.
    Other,
}

impl ObservedAllele {
    /// The allele a read's base stands for. Case-insensitive; everything else is
    /// [`Other`](Self::Other).
    pub fn of_base(base: u8) -> Self {
        match base.to_ascii_uppercase() {
            b'A' => Self::A,
            b'C' => Self::C,
            b'G' => Self::G,
            b'T' => Self::T,
            _ => Self::Other,
        }
    }

    pub fn code(self) -> u8 {
        match self {
            Self::A => 0,
            Self::C => 1,
            Self::G => 2,
            Self::T => 3,
            Self::Other => 4,
        }
    }
}

/// A depth, or the one state a depth cannot express.
///
/// A bin alone cannot say *this position was never visited*, and that has to be
/// distinguishable from *visited and empty* because only the first is a bug. So the stored
/// code is the ladder's bins plus one sentinel — 21 codes, which is why five bits and not
/// four.
///
/// **The ladder is not defined here.** It is [`DepthBinEdges`], shared with the histogram
/// route so the two cannot bin differently — which is what makes the comparison between the
/// routes a comparison.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum DepthCode {
    /// The region holding this position was never walked. **A bug, not data.**
    NeverWalked,
    /// A rung of the shared ladder; bin 0 is zero depth.
    Binned(DepthBin),
}

/// How many codes fit in one entry. Twenty bins plus the sentinel is 21, so five bits.
pub const DEPTH_CODE_BITS: u32 = 5;

/// The value [`DepthCode::NeverWalked`] is stored as — the top of the five-bit range, so
/// adding a rung to the ladder collides with it loudly rather than shifting it.
const NEVER_WALKED_CODE: u8 = (1 << DEPTH_CODE_BITS) - 1;

impl DepthCode {
    fn to_bits(self) -> u8 {
        match self {
            Self::NeverWalked => NEVER_WALKED_CODE,
            Self::Binned(bin) => {
                let raw = bin.get();
                assert!(
                    u32::from(raw) < NEVER_WALKED_CODE as u32,
                    "depth bin {raw} does not fit in {DEPTH_CODE_BITS} bits beside the \
                     never-walked sentinel; the ladder has outgrown the encoding"
                );
                raw as u8
            }
        }
    }

    fn from_bits(bits: u8) -> Self {
        if bits == NEVER_WALKED_CODE {
            Self::NeverWalked
        } else {
            Self::Binned(DepthBin(u16::from(bits)))
        }
    }
}

/// One depth code per kept position, five bits each, in selection order.
///
/// **No coordinates and no index**: the positions are reproducible from the selection rule,
/// so entry `i` is the `i`-th kept position and nothing says so on the wire. Storing
/// coordinates would cost about five bytes each — 10 MB at two million positions, before any
/// data was recorded.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PackedDepthCodes {
    bits: Vec<u8>,
    len: usize,
}

impl PackedDepthCodes {
    /// `len` entries, every one [`DepthCode::NeverWalked`].
    ///
    /// **That is the right initial state, not an empty one.** A kept position with no entry
    /// and a kept position never visited are the same thing, and the sentinel is what says
    /// so after the walk has finished and cannot be asked again.
    pub fn never_walked(len: usize) -> Self {
        let bytes = len.saturating_mul(DEPTH_CODE_BITS as usize).div_ceil(8);
        let mut packed = Self {
            bits: vec![0; bytes],
            len,
        };
        for index in 0..len {
            packed.set(index, DepthCode::NeverWalked);
        }
        packed
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// # Panics
    ///
    /// When `index` is past the end — a record set is sized from the selection, so an index
    /// outside it means the writer and the selection disagree, and that is not something to
    /// absorb.
    pub fn set(&mut self, index: usize, code: DepthCode) {
        assert!(
            index < self.len,
            "position {index} is outside the selection"
        );
        let (byte, offset) = (
            index * DEPTH_CODE_BITS as usize / 8,
            index * DEPTH_CODE_BITS as usize % 8,
        );
        let value = u32::from(code.to_bits());
        let mask = u32::from(NEVER_WALKED_CODE) << offset;
        let mut window = u32::from(self.bits[byte]);
        if byte + 1 < self.bits.len() {
            window |= u32::from(self.bits[byte + 1]) << 8;
        }
        window = (window & !mask) | (value << offset);
        self.bits[byte] = (window & 0xff) as u8;
        if byte + 1 < self.bits.len() {
            self.bits[byte + 1] = ((window >> 8) & 0xff) as u8;
        }
    }

    pub fn get(&self, index: usize) -> DepthCode {
        assert!(
            index < self.len,
            "position {index} is outside the selection"
        );
        let (byte, offset) = (
            index * DEPTH_CODE_BITS as usize / 8,
            index * DEPTH_CODE_BITS as usize % 8,
        );
        let mut window = u32::from(self.bits[byte]);
        if byte + 1 < self.bits.len() {
            window |= u32::from(self.bits[byte + 1]) << 8;
        }
        DepthCode::from_bits(((window >> offset) & u32::from(NEVER_WALKED_CODE)) as u8)
    }

    /// The packed bytes — 1.25 MB at two million positions.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bits
    }

    pub fn iter(&self) -> impl Iterator<Item = DepthCode> + '_ {
        (0..self.len).map(|index| self.get(index))
    }
}

/// One position's reads on one non-reference allele.
///
/// **Sparse, because at three reads a site nearly every position is "n reads, all matching"**
/// and what fills this list is the *error rate* rather than the variants: 30–250 kB at two
/// million positions.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct AlleleObservation {
    /// Index into the generic selection — the position's only identity.
    pub index: u32,
    pub allele: ObservedAllele,
    pub reads: u32,
}

/// One read group's generic evidence: a dense array of depths and a sparse list of what was
/// not on the reference base.
///
/// **The dense half reconstructs a quiet position exactly**: a code with no sparse entry
/// means that depth with every read on the reference base, and the reference base is a
/// property of the position rather than of the sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GenericEvidence {
    depth: PackedDepthCodes,
    non_reference: Vec<AlleleObservation>,
}

impl GenericEvidence {
    /// Records assembled from the two halves directly — **the door a reader comes in
    /// through, and the one a test that draws its own evidence uses.**
    ///
    /// # Panics
    ///
    /// When the sparse entries are not sorted by position, or name a position the depth array
    /// does not have. Both would leave the fit reading one position's reads at another's, and
    /// neither has a symptom.
    pub fn from_parts(depth: PackedDepthCodes, non_reference: Vec<AlleleObservation>) -> Self {
        assert!(
            non_reference
                .windows(2)
                .all(|pair| pair[0].index <= pair[1].index),
            "the sparse entries arrive in position order, and the fit walks them with a cursor"
        );
        assert!(
            non_reference
                .last()
                .is_none_or(|entry| (entry.index as usize) < depth.len()),
            "a sparse entry names a position the depth array does not have"
        );
        Self {
            depth,
            non_reference,
        }
    }

    pub fn never_walked(positions: usize) -> Self {
        Self {
            depth: PackedDepthCodes::never_walked(positions),
            non_reference: Vec::new(),
        }
    }

    pub fn depth(&self) -> &PackedDepthCodes {
        &self.depth
    }

    /// Ascending by `index`, then by allele — so the fit walks the dense and the sparse
    /// halves in one pass.
    pub fn non_reference(&self) -> &[AlleleObservation] {
        &self.non_reference
    }

    /// What this read group saw at the `index`-th kept position, as `(depth, alleles)`.
    pub fn at(&self, index: usize) -> (DepthCode, Vec<AlleleObservation>) {
        let start = self
            .non_reference
            .partition_point(|entry| (entry.index as usize) < index);
        let alleles = self.non_reference[start..]
            .iter()
            .take_while(|entry| entry.index as usize == index)
            .copied()
            .collect();
        (self.depth.get(index), alleles)
    }
}

// ---------------------------------------------------------------------
// The STR half of the census
// ---------------------------------------------------------------------

/// How far either side of the reference tract length a read's offset is recorded.
///
/// **Narrow is fine only because the end buckets are scored by their marginal**: "at least
/// four repeats short" gets the sum over every offset it absorbs, never the probability of
/// sitting exactly on the edge.
///
/// **Not the same constant as the span the fit may place allele mass on**, which reaches ±6.
/// That is the load-bearing one — it is what lets an end bucket be attributed to a distant
/// allele rather than to a far slip — and reading this width as that span is the confusion
/// the two constants exist to prevent.
pub const RECORDED_OFFSET_RANGE: i32 = 4;

/// Buckets, `-4 … +4` with saturating ends.
pub const OFFSET_BUCKETS: usize = (2 * RECORDED_OFFSET_RANGE + 1) as usize;

/// Spanning reads at each whole-repeat offset from the reference tract length.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct OffsetCounts {
    counts: [u16; OFFSET_BUCKETS],
}

impl Default for OffsetCounts {
    fn default() -> Self {
        Self {
            counts: [0; OFFSET_BUCKETS],
        }
    }
}

impl OffsetCounts {
    /// Add `reads` at `offset` whole repeat units from the reference tract length.
    ///
    /// **An offset beyond the recorded range saturates into the end bucket** rather than
    /// being dropped or wrapping — silent either way, and the end bucket is the only one the
    /// fit can score by a marginal.
    pub fn add(&mut self, offset: i32, reads: u16) {
        let clamped = offset.clamp(-RECORDED_OFFSET_RANGE, RECORDED_OFFSET_RANGE);
        let slot = (clamped + RECORDED_OFFSET_RANGE) as usize;
        self.counts[slot] = self.counts[slot].saturating_add(reads);
    }

    pub fn at(&self, offset: i32) -> u16 {
        let clamped = offset.clamp(-RECORDED_OFFSET_RANGE, RECORDED_OFFSET_RANGE);
        self.counts[(clamped + RECORDED_OFFSET_RANGE) as usize]
    }

    /// Reads that crossed the tract, whatever length they reported.
    pub fn total(&self) -> u32 {
        self.counts.iter().map(|c| u32::from(*c)).sum()
    }
}

/// A read whose tract differs from the reference by something that is not a whole number of
/// motif copies.
///
/// **A diagnostic, not a parameter.** Such a read is modelled as an independent per-read
/// outcome, so the likelihood splits exactly into how many reads were non-whole-repeat times
/// how the rest fell across the offsets; nothing about the slippage parameters is estimated
/// from it. **What it caught is recorded and not only how many**, because a partial unit at a
/// tract edge is alignment ambiguity and an indel in the flank is a different thing, and a
/// bare count can raise a threshold without ever explaining it.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct GuardObservation {
    /// Index into the STR selection.
    pub locus: u32,
    /// The read's tract length minus the reference's, in bases.
    pub length_difference: i32,
    pub reads: u16,
}

/// One mismatching base on one read.
///
/// **Offsets record length, and a substitution that does not change a tract's length is
/// invisible to them**, so without this channel the STR error rate cannot be recovered at
/// all. A pair of counters would do for the rate and for nothing else: it cannot separate a
/// substitution inside the tract — which interrupts the motif and changes what a repeat unit
/// is — from ordinary error in the flank, and it cannot say that two reads carried the *same*
/// interruption, which is what makes an interruption an allele rather than an error.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct TractDifference {
    /// Index into the STR selection.
    pub locus: u32,
    /// Which of this locus's reads carried it, in the locus's own read order.
    /// **Two entries at one offset on one read is a different observation from the same two
    /// on two reads**, and a read-blind encoding passes every other check in the suite.
    pub read: u8,
    /// Signed offset from the tract start: negative in the left flank, `0..len` inside the
    /// tract, at or beyond `len` in the right flank.
    pub offset: i16,
    pub base: ObservedAllele,
}

/// One read group's STR evidence.
///
/// **Four states at a locus, not two**: no read reached it; reads reached it but none
/// crossed the whole tract, so none reports a length; reads crossed it; and the region was
/// never walked. The generic set's zero-depth-against-quiet distinction is the last pair, and
/// [`covering_not_crossing`](Self::covering_not_crossing) is what makes the second
/// expressible.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrEvidence {
    offsets: Vec<OffsetCounts>,
    covering_not_crossing: Vec<u16>,
    walked: Vec<bool>,
    bases_compared: Vec<u32>,
    guard: Vec<GuardObservation>,
    differences: Vec<TractDifference>,
}

impl SsrEvidence {
    pub fn never_walked(loci: usize) -> Self {
        Self {
            offsets: vec![OffsetCounts::default(); loci],
            covering_not_crossing: vec![0; loci],
            walked: vec![false; loci],
            bases_compared: vec![0; loci],
            guard: Vec::new(),
            differences: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    pub fn offsets(&self, locus: usize) -> OffsetCounts {
        self.offsets[locus]
    }

    /// Reads that reached this locus and crossed no whole tract — a censored lower bound on
    /// the length.
    ///
    /// **The censoring is not random**, which is why it is a field and not an inference: a
    /// tract longer than a read is never crossed, in every sample at every depth, so it runs
    /// along repeat count — the very axis the slippage numbers are fitted within — and a
    /// stratum unobservable with this read length must not look like one that was merely
    /// unlucky with coverage.
    pub fn covering_not_crossing(&self, locus: usize) -> u16 {
        self.covering_not_crossing[locus]
    }

    /// The denominator the STR error rate is fitted against.
    pub fn bases_compared(&self, locus: usize) -> u32 {
        self.bases_compared[locus]
    }

    pub fn guard(&self) -> &[GuardObservation] {
        &self.guard
    }

    pub fn differences(&self) -> &[TractDifference] {
        &self.differences
    }

    /// Which of the four states this locus is in.
    pub fn state(&self, locus: usize) -> SsrLocusState {
        if !self.walked[locus] {
            SsrLocusState::NeverWalked
        } else if self.offsets[locus].total() > 0 {
            SsrLocusState::Crossed
        } else if self.covering_not_crossing[locus] > 0 {
            SsrLocusState::ReachedNotCrossed
        } else {
            SsrLocusState::NoRead
        }
    }

    /// Whether the guard has caught more than one read in ten of those that differ from the
    /// reference tract length.
    ///
    /// **Above that the locus is not something this noise model describes**, and the fit
    /// should say so rather than fit it. It is not an error at write time: a locus over the
    /// threshold is well-formed data, and encoding it as a write failure would make a
    /// property of the sample look like a property of the file.
    pub fn guard_is_over_threshold(&self, locus: usize) -> bool {
        let guarded: u32 = self
            .guard
            .iter()
            .filter(|entry| entry.locus as usize == locus)
            .map(|entry| u32::from(entry.reads))
            .sum();
        let differing: u32 = (-RECORDED_OFFSET_RANGE..=RECORDED_OFFSET_RANGE)
            .filter(|offset| *offset != 0)
            .map(|offset| u32::from(self.offsets[locus].at(offset)))
            .sum();
        guarded * 10 > differing + guarded
    }
}

/// Which of the four states a sample is in at a kept STR locus.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SsrLocusState {
    /// The region was never walked — a bug, not data.
    NeverWalked,
    /// Walked, and no read reached the locus.
    NoRead,
    /// Reads reached the locus and none crossed the whole tract: a lower bound on the
    /// length, and no length reported.
    ReachedNotCrossed,
    /// At least one read crossed the tract.
    Crossed,
}

// ---------------------------------------------------------------------
// What travels beside the evidence
// ---------------------------------------------------------------------

/// How many of a locus's reads are entered, before the rest are subsampled away.
///
/// **A newtype and not a `usize`** because it travels in [`RecordingTerms`] beside other
/// counts and is compared for equality: a bare integer there is transposable with
/// [`DepthCap`], and the two are not the same number.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct ReadCap(pub u32);

/// The depth above which a position's reads are subsampled before anything is recorded.
///
/// The generic path's twin of [`ReadCap`]. **It moves independently of the ladder's own top
/// rung**, and a sample recorded at a different one did not record the same evidence.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct DepthCap(pub u32);

/// A digest of the depth ladder's edges.
///
/// **The generic record stores a five-bit code, not a depth.** Two samples binned under
/// different edges hold codes that mean different depths, *and every other value in
/// [`RecordingTerms`] agrees* — the loci were the same, the seed was the same, the digest of
/// the kept loci matches because the loci did match. A code is only a number until something
/// says what ladder it indexes.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct DepthLadderDigest(pub [u8; 16]);

impl DepthLadderDigest {
    pub fn of(edges: &DepthBinEdges) -> Self {
        let mut hasher = Md5::new();
        for bin in edges.bins() {
            let range = edges.depth_range(bin);
            hasher.update(range.start().to_le_bytes());
            hasher.update(range.end().to_le_bytes());
        }
        Self(hasher.finalize().into())
    }
}

/// The thirteen values the fit refuses to pool across.
///
/// Seven say which loci were **asked for** ([`SelectionTerms`]), one says which came
/// **back** ([`CensusLociDigest`]), and five say **in what units** the evidence was written
/// down. An earlier version of this check had none of the last five, which left it able to
/// pass while two samples' rows meant different things.
#[derive(Debug, Clone, PartialEq)]
pub struct RecordingTerms {
    pub selection: SelectionTerms,
    pub kept_loci: CensusLociDigest,
    /// Per stratum, how many loci the analysed regions hold against how many were kept.
    /// Anything pooled across strata is biased without it and silently so.
    pub ssr_stratum_counts: StratumCounts,
    pub read_cap: ReadCap,
    pub depth_ladder: DepthLadderDigest,
    pub depth_cap: DepthCap,
    /// The coverage-by-window grid, where a summary exists. Windows of different widths are
    /// not comparable and a relative copy number computed across two grids is meaningless.
    pub coverage_window: Option<Bp>,
}

impl RecordingTerms {
    /// Which value two samples first disagree on — `None` when they may be pooled.
    ///
    /// **Naming the value is the whole point.** Every value here fails the same way,
    /// silently, and only the name says what to fix; the selection's own values delegate to
    /// their own check for the same reason.
    ///
    /// **Destructured without `..` on purpose.** A value added to this struct stops this
    /// function compiling rather than quietly going unchecked, and a value that goes
    /// unchecked lets two samples that recorded different things be pooled without a word.
    pub fn first_disagreement(&self, other: &Self) -> Option<&'static str> {
        let Self {
            selection,
            kept_loci,
            ssr_stratum_counts,
            read_cap,
            depth_ladder,
            depth_cap,
            coverage_window,
        } = self;
        if let Some(field) = selection.first_disagreement(&other.selection) {
            return Some(field);
        }
        if kept_loci != &other.kept_loci {
            return Some("the loci actually kept");
        }
        if ssr_stratum_counts.iter_sorted() != other.ssr_stratum_counts.iter_sorted() {
            return Some("per-stratum locus counts");
        }
        if read_cap != &other.read_cap {
            return Some("per-locus read cap");
        }
        if depth_ladder != &other.depth_ladder {
            return Some("depth ladder edges");
        }
        if depth_cap != &other.depth_cap {
            return Some("per-position depth cap");
        }
        if coverage_window != &other.coverage_window {
            return Some("coverage window size");
        }
        None
    }
}

// ---------------------------------------------------------------------
// The third object, which is not census evidence
// ---------------------------------------------------------------------

/// One sample's depth over fixed windows of the reference, plus the GC curve that corrects
/// it. **Over every position the walk visited, not over the kept ones.**
///
/// **Why it cannot come out of the records.** The fit's third class of site — a locus the
/// sample carries more copies of than the reference does — is conditioned on local relative
/// coverage rather than on the site's own depth, and the records hold one binned depth per
/// kept position where the kept positions are one in a few hundred. A 500 bp window holds one
/// or two of them, which is the per-base measurement the class's own constraint rules out.
///
/// **Measured, 2026-08-12** (`doc/devel/ng/reports/duplicated_locus_probe_2026-08-12.md`): on
/// tomato at 25× depth, 1 position in 8,600 sits in a window near two copies and reads
/// between 35% and 65% alternative, and the near-half rate inside those windows is 1.26%
/// against 0.033% outside.
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageByWindow {
    window_bp: Bp,
    median_depth: f32,
    /// `round(32 × mean / median_depth)`, saturating at 255 — 3% of the sample's own median
    /// per step, reaching eight times it. **The mean itself in a byte would not do**: at
    /// three reads a position the difference between one copy and two is the difference
    /// between 3 and 6.
    depth: Vec<u8>,
    /// Windows holding fewer than `window_bp` walked positions, as `(index, positions)` —
    /// contig ends, and anything the analysed regions or the ambiguity mask cut into.
    short_windows: Vec<(u32, u16)>,
    gc_curve: Vec<f32>,
}

/// The scale a window's mean depth is stored on, relative to the sample's median.
const WINDOW_DEPTH_SCALE: f32 = 32.0;

impl CoverageByWindow {
    /// Build from one mean depth per window and the count of positions behind each.
    ///
    /// # Panics
    ///
    /// When the two slices differ in length — they are two descriptions of one grid.
    pub fn new(
        window_bp: Bp,
        median_depth: f32,
        means: &[f32],
        positions: &[u16],
        gc_curve: Vec<f32>,
    ) -> Self {
        assert_eq!(
            means.len(),
            positions.len(),
            "every window's mean needs the count of positions behind it"
        );
        let full = u16::try_from(window_bp.get()).unwrap_or(u16::MAX);
        let depth = means
            .iter()
            .map(|mean| {
                if median_depth <= 0.0 {
                    return 0;
                }
                let scaled = (WINDOW_DEPTH_SCALE * mean / median_depth).round();
                scaled.clamp(0.0, 255.0) as u8
            })
            .collect();
        let short_windows = positions
            .iter()
            .enumerate()
            .filter(|(_, count)| **count < full)
            .map(|(index, count)| (index as u32, *count))
            .collect();
        Self {
            window_bp,
            median_depth,
            depth,
            short_windows,
            gc_curve,
        }
    }

    pub fn window_bp(&self) -> Bp {
        self.window_bp
    }

    pub fn len(&self) -> usize {
        self.depth.len()
    }

    pub fn is_empty(&self) -> bool {
        self.depth.is_empty()
    }

    pub fn gc_curve(&self) -> &[f32] {
        &self.gc_curve
    }

    /// **The sample's own median window depth — the scale every relative reading is against.**
    ///
    /// A duplicated stretch is *twice this sample's normal*, never an absolute depth, so
    /// anything asking whether a window holds one copy or two divides by this first.
    #[must_use]
    pub fn median_depth(&self) -> f32 {
        self.median_depth
    }

    /// This window's mean depth, in reads a position.
    pub fn mean_depth(&self, window: usize) -> f32 {
        f32::from(self.depth[window]) * self.median_depth / WINDOW_DEPTH_SCALE
    }

    /// How many walked positions stand behind this window.
    pub fn positions(&self, window: usize) -> u16 {
        let full = u16::try_from(self.window_bp.get()).unwrap_or(u16::MAX);
        self.short_windows
            .binary_search_by_key(&(window as u32), |(index, _)| *index)
            .map_or(full, |slot| self.short_windows[slot].1)
    }

    /// The mean depth over `windows` taken together — **the operation the fit reads through**.
    ///
    /// **The stored grid is fine and the width read at is the sample's own.** A window's mean
    /// separates one copy from two only once it has collected about 12,000 aligned bases, so
    /// 500 bp is enough at 25× and 5 kb is needed at 2.5×; at 3.6× and 500 bp the enrichment
    /// of the joint cell over independence falls to 1.3, which is no separation at all
    /// (`parameter_prepass_joint_records.md` §4.1).
    ///
    /// It is a ratio of two sums, so it is exact — and it weights each window by its own
    /// position count, which is why those counts are stored. A sum that treated every window
    /// as full would agree everywhere except at contig ends and region edges, which is where
    /// it matters.
    pub fn mean_depth_over(&self, windows: std::ops::Range<usize>) -> f32 {
        let mut depth_sum = 0.0_f64;
        let mut positions = 0.0_f64;
        for window in windows {
            let count = f64::from(self.positions(window));
            depth_sum += f64::from(self.mean_depth(window)) * count;
            positions += count;
        }
        if positions == 0.0 {
            0.0
        } else {
            (depth_sum / positions) as f32
        }
    }

    /// How many adjacent windows this sample's depth needs summed before a window's mean
    /// separates one copy from two.
    ///
    /// [`MIN_ALIGNED_BASES_PER_WINDOW`] over the sample's median depth and the stored width,
    /// rounded up, and never below one.
    pub fn windows_to_sum(&self) -> usize {
        if self.median_depth <= 0.0 {
            return 1;
        }
        let per_window = f64::from(self.median_depth) * self.window_bp.get() as f64;
        ((f64::from(MIN_ALIGNED_BASES_PER_WINDOW) / per_window).ceil() as usize).max(1)
    }
}

/// How many aligned bases a window has to collect before its mean depth tells one copy from
/// two.
///
/// **Measured across eight tomato samples from 2.5× to 28.7×**: 2.51× at 5 kb is 12,550
/// aligned bases a window and 25.2× at 500 bp is 12,600, and the two return the same
/// enrichment of the joint cell — 14× and 24× — while 3.6× at 500 bp, which is 1,800, returns
/// 1.3 and separates nothing. The deep sample gains nothing above the floor, so this is a
/// floor and not a target.
pub const MIN_ALIGNED_BASES_PER_WINDOW: u32 = 12_000;

// ---------------------------------------------------------------------
// The whole input to the fit
// ---------------------------------------------------------------------

/// One sample's evidence at the kept loci, plus the values the fit checks before pooling.
///
/// **`BTreeMap` and not `HashMap`**, so a fit that iterates read groups is deterministic. A
/// sample's counts at a position are the sum of its read groups' — raw counts at one place,
/// so the equality is exact and the fit may fold freely.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleCensusEvidence {
    pub sample: String,
    pub generic: BTreeMap<ReadGroupId, GenericEvidence>,
    pub ssr: BTreeMap<ReadGroupId, SsrEvidence>,
    pub coverage: Option<CoverageByWindow>,
    pub terms: RecordingTerms,
}

impl SampleCensusEvidence {
    /// This sample's depth at the `index`-th kept generic position, summed over its read
    /// groups.
    ///
    /// **Summing is where the two grains meet.** The cross-sample work wants the *sample's*
    /// counts at a position; the error rate is chemistry and wants the *read group's*. The
    /// record is kept at the finer grain because summing is addition of raw counts at one
    /// place and the reverse is not recoverable.
    pub fn allele_counts(&self, index: usize) -> [u32; 5] {
        let mut counts = [0_u32; 5];
        for records in self.generic.values() {
            for observation in records.at(index).1 {
                counts[observation.allele.code() as usize] += observation.reads;
            }
        }
        counts
    }
}

// ---------------------------------------------------------------------
// Filling the census during the walk
// ---------------------------------------------------------------------

/// Fills a sample's records as the walk hands it loci.
///
/// **It is handed the same locus stream the histogram accumulators are handed**, borrowing
/// rather than taking, so one walk fills both routes and the comparison between them is over
/// identical evidence. It knows which loci are kept and ignores the rest.
pub struct CensusWriter {
    /// The kept generic positions, in selection order — the index a record entry carries.
    generic_loci: Vec<GenomePosition>,
    /// The kept STR loci as regions, in genome order.
    ssr_loci: Vec<GenomeRegion>,
    edges: DepthBinEdges,
    generic: BTreeMap<ReadGroupId, GenericEvidence>,
    ssr: BTreeMap<ReadGroupId, SsrEvidence>,
    /// Non-reference observations before they are sorted — the walk arrives in position
    /// order, so this is already sorted in practice and the sort at `finish` is a guard.
    pending: BTreeMap<ReadGroupId, Vec<AlleleObservation>>,
    digester: CensusLociDigester,
    /// How far the digester has been fed, so a locus arriving out of order is caught.
    digested: usize,
    terms: SelectionTerms,
    read_cap: ReadCap,
    depth_cap: DepthCap,
    stratum_counts: StratumCounts,
    /// Every read group the sample declares. **Held, rather than discovered from the
    /// observations**, so a group that put no read at a position still gets its zero there:
    /// discovering them would make a group's record start at its first read, and every
    /// position before that indistinguishable from never walked.
    read_groups: Vec<ReadGroupId>,
    sample: String,
}

/// How many of `reads` survive when a position's `depth` reads are thinned down to `cap`.
///
/// **Proportional and rounded to nearest, never a draw.** A region-sharded walk and a
/// single-threaded one must produce byte-identical records, and a share needs no seed to do
/// that. Rounding to nearest keeps a single stray read at a very deep position rather than
/// discarding it, which matters because those reads are what the error rate is fitted from.
fn thin_to_cap(reads: u32, depth: u32, cap: u32) -> u32 {
    if depth <= cap || depth == 0 {
        return reads;
    }
    let scaled = f64::from(reads) * f64::from(cap) / f64::from(depth);
    (scaled.round() as u32).min(cap).max(u32::from(reads > 0))
}

impl CensusWriter {
    #[allow(
        clippy::too_many_arguments,
        reason = "every one is a distinct input the records cannot be built without; \
                  bundling them into a config struct would only move the list"
    )]
    pub fn new(
        sample: String,
        loci: &CensusLoci,
        read_groups: Vec<ReadGroupId>,
        contig_of: &dyn Fn(&str) -> Option<ContigId>,
        terms: SelectionTerms,
        edges: DepthBinEdges,
        read_cap: ReadCap,
        depth_cap: DepthCap,
    ) -> Self {
        let mut ssr_loci: Vec<GenomeRegion> = loci
            .ssr()
            .iter_sorted()
            .into_iter()
            .flat_map(|(_, segments)| segments.iter())
            .filter_map(|segment| {
                contig_of(segment.chrom()).map(|contig| GenomeRegion {
                    contig,
                    start: Position(segment.start()),
                    end: Position(segment.end()),
                })
            })
            .collect();
        ssr_loci.sort_unstable_by_key(|r| (r.contig.get(), r.start.get(), r.end.get()));
        Self {
            generic_loci: loci.generic().to_vec(),
            ssr: BTreeMap::new(),
            ssr_loci,
            edges,
            generic: BTreeMap::new(),
            pending: BTreeMap::new(),
            digester: CensusLociDigester::new(),
            digested: 0,
            terms,
            read_cap,
            depth_cap,
            stratum_counts: loci.ssr_stratum_counts().clone(),
            read_groups,
            sample,
        }
    }

    /// Record this locus if it is a kept one.
    ///
    /// **Every kept locus gets an entry whether or not a read reached it** — the entry is the
    /// denominator, and a locus in a region never walked keeps [`DepthCode::NeverWalked`], so
    /// the three states survive.
    pub fn add_locus(&mut self, locus: &SampleLocusObservations) {
        match &locus.kind {
            LocusKind::Generic => self.add_generic(locus),
            LocusKind::Ssr(_) => self.add_ssr(locus),
            LocusKind::SsrBundle => {}
        }
    }

    /// This stretch of genome was walked, whatever the walk found in it.
    ///
    /// **Without this the three states of §1.1 collapse to two on real data.** The generic
    /// locus generator emits a locus only where a read reached, so a position no read reached
    /// produces nothing and would keep [`DepthCode::NeverWalked`] — which is the code for a
    /// bug, a region the run never opened. Measured on tomato SRR7279482 at 25×, that is
    /// 93,150 of 1,999,404 kept positions, 1 in 21, every one of them data being reported as
    /// a defect.
    ///
    /// Call it for each region handed to the walk, before or after the loci from it: a real
    /// depth overwrites a zero and a zero never overwrites a real depth.
    pub fn mark_walked(&mut self, region: GenomeRegion) {
        let zero = DepthCode::Binned(self.edges.bin_for(0));
        let first = self.generic_loci.partition_point(|kept| {
            (kept.contig.get(), kept.position.get()) < (region.contig.get(), region.start.get())
        });
        for index in first..self.generic_loci.len() {
            let kept = self.generic_loci[index];
            if kept.contig != region.contig || kept.position.get() > region.end.get() {
                break;
            }
            for group in &self.read_groups {
                let records = self
                    .generic
                    .entry(*group)
                    .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
                if records.depth.get(index) == DepthCode::NeverWalked {
                    records.depth.set(index, zero);
                }
            }
        }
    }

    fn add_generic(&mut self, locus: &SampleLocusObservations) {
        // **Depth is per read group, not pooled.** `num_obs_along_locus` sums the sample's
        // observations, which is the wrong grain for a record keyed by read group — and it
        // is the grain the error rate is fitted at, so pooling here would score every read
        // against a rate fitted from a depth its own library never had.
        let span = locus.region.len() as usize;
        let mut depths: BTreeMap<ReadGroupId, Vec<u32>> = BTreeMap::new();
        for group in &self.read_groups {
            depths.insert(*group, vec![0; span]);
        }
        for observation in &locus.observations {
            let Some(depth) = depths.get_mut(&observation.read_group) else {
                continue;
            };
            match &observation.read_witness {
                ReadWitness::Complete => {
                    for slot in depth.iter_mut() {
                        *slot = slot.saturating_add(observation.num_obs);
                    }
                }
                ReadWitness::Partial { positions } => {
                    for (start, end) in positions.runs() {
                        let from = (start as usize).min(span);
                        let to = (end as usize).min(span).max(from);
                        for slot in &mut depth[from..to] {
                            *slot = slot.saturating_add(observation.num_obs);
                        }
                    }
                }
            }
        }

        // The whole locus's observations describe the whole locus's span; only a
        // single-base locus lets "which base did this read show here" be answered from them.
        // A wider one contributes its depth, which is what the fit reads it for.
        let single_base = span == 1;
        for offset in 0..span {
            let position = GenomePosition {
                contig: locus.region.contig,
                position: Position(locus.region.start.get() + offset as u64),
            };
            let Some(index) = self.generic_index(position) else {
                continue;
            };
            self.feed_digest(index, position);

            // **Every read group gets an entry at every walked position**, whether or not it
            // put a read there: the entry is the denominator, and a group with no read here
            // saw zero, which is data and not the absence of it.
            //
            // **The depth cap is applied here and nowhere else** (§5). Above it the reads are
            // thinned to it, keeping the fractions they showed, because the stored code is a
            // bin and the ladder's top bin is the last one there is: a sample at 300 reads a
            // position whose reads were *not* thinned records a depth of about 111 beside an
            // undiminished count of alternative reads, which charges it a negative number of
            // reference reads and inflates every rate fitted from it.
            for (group, depth) in &depths {
                let records = self
                    .generic
                    .entry(*group)
                    .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
                let capped = depth[offset].min(self.depth_cap.0);
                records
                    .depth
                    .set(index, DepthCode::Binned(self.edges.bin_for(capped)));
            }
            if !single_base {
                continue;
            }
            for observation in locus.complete_observations() {
                if *observation.bases == *locus.reference_bases {
                    continue;
                }
                let allele = match &*observation.bases {
                    [base] => ObservedAllele::of_base(*base),
                    // An insertion or a deletion at a one-base locus: not one of the four,
                    // and the fit scores it as the fifth rather than guessing a base.
                    _ => ObservedAllele::Other,
                };
                // **Thinned by the same ratio as the depth**, so the fractions the reads showed
                // survive the cap exactly. Deterministic rather than a draw: a region-sharded
                // walk and a single-threaded one must produce byte-identical records, and a
                // proportional share needs no seed to do that.
                let group_depth = depths
                    .get(&observation.read_group)
                    .map_or(0, |depth| depth[offset]);
                let reads = thin_to_cap(observation.num_obs, group_depth, self.depth_cap.0);
                if reads == 0 {
                    continue;
                }
                self.pending
                    .entry(observation.read_group)
                    .or_default()
                    .push(AlleleObservation {
                        index: index as u32,
                        allele,
                        reads,
                    });
            }
        }
    }

    fn add_ssr(&mut self, locus: &SampleLocusObservations) {
        let Some(index) = self.ssr_index(locus.region) else {
            return;
        };
        let LocusKind::Ssr(detail) = &locus.kind else {
            return;
        };
        let period = detail.motif.period().max(1) as i64;
        let reference_length = locus.region.len() as i64;

        for observation in &locus.observations {
            let records = self
                .ssr
                .entry(observation.read_group)
                .or_insert_with(|| SsrEvidence::never_walked(self.ssr_loci.len()));
            records.walked[index] = true;
            let reads = u16::try_from(observation.num_obs).unwrap_or(u16::MAX);
            if observation.read_witness != crate::ng::locus_generation::ReadWitness::Complete {
                // A read that covered the tract without crossing it reports no length; it is
                // a lower bound, and the fit is told so rather than shown a short allele.
                records.covering_not_crossing[index] =
                    records.covering_not_crossing[index].saturating_add(reads);
                continue;
            }
            let difference = observation.bases.len() as i64 - reference_length;
            if difference % period == 0 {
                records.offsets[index].add((difference / period) as i32, reads);
            } else {
                records.guard.push(GuardObservation {
                    locus: index as u32,
                    length_difference: difference as i32,
                    reads,
                });
            }
            // **The difference list, and the one case it can be built from.** A read whose
            // tract is the reference's length lines up with it base for base, so a mismatch is
            // read off directly. A read that slipped does not: which of its bases sits over
            // which reference base is the aligner's answer and not this writer's, and inventing
            // one would put a whole tract's worth of manufactured mismatches into the error
            // rate. Such a read contributes its length to the offsets and **nothing to the
            // denominator**, so the rate stays a ratio of two quantities counted over the same
            // reads.
            if observation.bases.len() as i64 != reference_length {
                continue;
            }
            records.bases_compared[index] +=
                u32::from(reads) * u32::try_from(observation.bases.len()).unwrap_or(u32::MAX);
            // **One entry per read, not one per distinct sequence.** The walk has already
            // folded reads carrying the same bases into a single observation with a count, and
            // collapsing them here too would lose the fact this channel exists for: that two
            // reads carried the *same* interruption, which is what makes an interruption an
            // allele rather than an error.
            for (offset, (read_base, reference_base)) in observation
                .bases
                .iter()
                .zip(locus.reference_bases.iter())
                .enumerate()
            {
                if read_base.eq_ignore_ascii_case(reference_base) {
                    continue;
                }
                for copy in 0..reads {
                    records.differences.push(TractDifference {
                        locus: index as u32,
                        read: u8::try_from(copy).unwrap_or(u8::MAX),
                        offset: i16::try_from(offset).unwrap_or(i16::MAX),
                        base: ObservedAllele::of_base(*read_base),
                    });
                }
            }
        }
        // Reads that reached the locus and produced no observation at all.
        if locus.reads_without_observation > 0 {
            for records in self.ssr.values_mut() {
                records.walked[index] = true;
            }
            if let Some(records) = self.ssr.values_mut().next() {
                let reads = u16::try_from(locus.reads_without_observation).unwrap_or(u16::MAX);
                records.covering_not_crossing[index] =
                    records.covering_not_crossing[index].saturating_add(reads);
            }
        }
    }

    fn feed_digest(&mut self, index: usize, position: GenomePosition) {
        // The walk arrives in genome order, so the digest is fed in index order; a locus
        // reached twice is fed once.
        while self.digested <= index {
            let at = self.digested;
            self.digester.observe(at, self.generic_loci[at]);
            self.digested += 1;
        }
        debug_assert_eq!(self.generic_loci[index], position);
    }

    fn generic_index(&self, position: GenomePosition) -> Option<usize> {
        self.generic_loci
            .binary_search_by_key(&(position.contig.get(), position.position.get()), |kept| {
                (kept.contig.get(), kept.position.get())
            })
            .ok()
    }

    fn ssr_index(&self, region: GenomeRegion) -> Option<usize> {
        self.ssr_loci
            .binary_search_by_key(&(region.contig.get(), region.start.get()), |kept| {
                (kept.contig.get(), kept.start.get())
            })
            .ok()
    }

    /// The finished records, with the digest of what was actually written.
    pub fn finish(mut self, coverage: Option<CoverageByWindow>) -> SampleCensusEvidence {
        // Every kept locus is digested, including any the walk never reached — the digest
        // witnesses the selection, and a run that stopped early must not produce a short
        // digest that happens to match another short one.
        while self.digested < self.generic_loci.len() {
            let at = self.digested;
            self.digester.observe(at, self.generic_loci[at]);
            self.digested += 1;
        }
        // **Every declared read group ends with a record set, even one that recorded
        // nothing.** A missing map entry and a set of never-walked entries are the same
        // situation, and only the second says so — a fit handed the first would have to
        // guess whether the group exists.
        for group in &self.read_groups {
            self.generic
                .entry(*group)
                .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
        }
        for (group, mut entries) in std::mem::take(&mut self.pending) {
            entries.sort_unstable_by_key(|entry| (entry.index, entry.allele));
            let records = self
                .generic
                .entry(group)
                .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
            records.non_reference = entries;
        }
        let coverage_window = coverage.as_ref().map(CoverageByWindow::window_bp);
        SampleCensusEvidence {
            sample: self.sample,
            generic: self.generic,
            ssr: self.ssr,
            coverage,
            terms: RecordingTerms {
                selection: self.terms,
                kept_loci: self.digester.finish(),
                ssr_stratum_counts: self.stratum_counts,
                read_cap: self.read_cap,
                depth_ladder: DepthLadderDigest::of(&self.edges),
                depth_cap: self.depth_cap,
                coverage_window,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- the depth code, and the state a depth cannot express ----------------

    #[test]
    fn every_code_survives_packing_at_every_bit_offset() {
        let edges = DepthBinEdges::new();
        let mut codes: Vec<DepthCode> = edges.bins().map(DepthCode::Binned).collect();
        codes.push(DepthCode::NeverWalked);
        // Repeat the ladder so that each code lands at every offset within a byte.
        let written: Vec<DepthCode> = (0..8).flat_map(|_| codes.clone()).collect();

        let mut packed = PackedDepthCodes::never_walked(written.len());
        for (index, code) in written.iter().enumerate() {
            packed.set(index, *code);
        }
        assert_eq!(packed.iter().collect::<Vec<_>>(), written);
    }

    #[test]
    fn a_fresh_array_is_never_walked_rather_than_zero_depth() {
        // The three states the spec requires, and the one this checks is the bug.
        let packed = PackedDepthCodes::never_walked(3);
        assert!(packed.iter().all(|code| code == DepthCode::NeverWalked));
        assert_ne!(
            DepthCode::NeverWalked,
            DepthCode::Binned(DepthBinEdges::new().bin_for(0))
        );
    }

    #[test]
    fn five_bits_a_position_and_no_more() {
        // 1.25 MB at two million positions is the size the spec prices, and it is this.
        let packed = PackedDepthCodes::never_walked(2_000_000);
        assert_eq!(packed.as_bytes().len(), 1_250_000);
    }

    // ---- the generic half ---------------------------------------------------

    fn generic_evidence() -> GenericEvidence {
        // Five positions, one per state: walked at zero depth, reads with a non-reference
        // allele, reads with two of them, reads with none at all, and never walked.
        let mut records = GenericEvidence::never_walked(5);
        let edges = DepthBinEdges::new();
        records.depth.set(0, DepthCode::Binned(edges.bin_for(0)));
        records.depth.set(1, DepthCode::Binned(edges.bin_for(7)));
        records.depth.set(2, DepthCode::Binned(edges.bin_for(300)));
        records.depth.set(3, DepthCode::Binned(edges.bin_for(5)));
        records.non_reference = vec![
            AlleleObservation {
                index: 1,
                allele: ObservedAllele::C,
                reads: 2,
            },
            AlleleObservation {
                index: 2,
                allele: ObservedAllele::G,
                reads: 1,
            },
            AlleleObservation {
                index: 2,
                allele: ObservedAllele::Other,
                reads: 3,
            },
        ];
        records
    }

    #[test]
    fn the_dense_half_reconstructs_a_quiet_position_and_the_sparse_half_the_rest() {
        let records = generic_evidence();
        let edges = DepthBinEdges::new();

        // Walked, zero depth: data, and distinguishable from never walked.
        assert_eq!(records.at(0), (DepthCode::Binned(edges.bin_for(0)), vec![]));
        // Reads, none of them non-reference: also data, and a third thing again.
        assert_eq!(
            records.at(3),
            (DepthCode::Binned(edges.bin_for(5)), vec![]),
            "five reads and no sparse entry means five reads on the reference base"
        );
        // Never walked, which is the bug the other two must not be confused with.
        assert_eq!(records.at(4).0, DepthCode::NeverWalked);
        // A multi-allelic position comes back with both alleles.
        let (depth, alleles) = records.at(2);
        assert_eq!(depth, DepthCode::Binned(edges.bin_for(300)));
        assert_eq!(alleles.len(), 2);
        assert_eq!(alleles[0].allele, ObservedAllele::G);
        assert_eq!(alleles[1].allele, ObservedAllele::Other);
    }

    #[test]
    fn read_groups_fold_by_addition() {
        // Two read groups, the same position: the sample's count is the sum, exactly.
        let mut one = GenericEvidence::never_walked(2);
        one.non_reference = vec![AlleleObservation {
            index: 1,
            allele: ObservedAllele::C,
            reads: 2,
        }];
        let mut two = GenericEvidence::never_walked(2);
        two.non_reference = vec![AlleleObservation {
            index: 1,
            allele: ObservedAllele::C,
            reads: 5,
        }];
        let records = SampleCensusEvidence {
            sample: "s".to_string(),
            generic: BTreeMap::from([(ReadGroupId(0), one), (ReadGroupId(1), two)]),
            ssr: BTreeMap::new(),
            coverage: None,
            terms: terms(),
        };
        assert_eq!(records.allele_counts(1), [0, 7, 0, 0, 0]);
        assert_eq!(
            records.allele_counts(0),
            [0; 5],
            "the fold reads the position it was asked for, not every entry there is"
        );
    }

    // ---- the STR half -------------------------------------------------------

    #[test]
    fn an_offset_past_the_recorded_range_saturates_rather_than_wrapping() {
        let mut counts = OffsetCounts::default();
        counts.add(7, 3);
        counts.add(-9, 2);
        counts.add(0, 1);
        assert_eq!(counts.at(RECORDED_OFFSET_RANGE), 3);
        assert_eq!(counts.at(-RECORDED_OFFSET_RANGE), 2);
        assert_eq!(counts.at(0), 1);
        assert_eq!(counts.total(), 6);
        // The end bucket holds the marginal: everything at or beyond it, not the edge alone.
        counts.add(4, 5);
        assert_eq!(counts.at(RECORDED_OFFSET_RANGE), 8);
    }

    #[test]
    fn the_four_states_at_an_str_locus_are_distinguishable() {
        let mut records = SsrEvidence::never_walked(4);
        assert_eq!(records.state(0), SsrLocusState::NeverWalked);

        records.walked[1] = true;
        assert_eq!(records.state(1), SsrLocusState::NoRead);

        records.walked[2] = true;
        records.covering_not_crossing[2] = 3;
        assert_eq!(records.state(2), SsrLocusState::ReachedNotCrossed);

        records.walked[3] = true;
        records.offsets[3].add(0, 5);
        assert_eq!(records.state(3), SsrLocusState::Crossed);
    }

    #[test]
    fn the_difference_list_tells_one_interruption_on_two_reads_from_two_errors() {
        // The same interior offset on two reads is an allele; two offsets on one read is
        // not. A read-blind encoding cannot tell them apart and passes every other check.
        let shared = vec![
            TractDifference {
                locus: 0,
                read: 0,
                offset: 3,
                base: ObservedAllele::A,
            },
            TractDifference {
                locus: 0,
                read: 1,
                offset: 3,
                base: ObservedAllele::A,
            },
        ];
        let scattered = vec![
            TractDifference {
                locus: 0,
                read: 0,
                offset: 3,
                base: ObservedAllele::A,
            },
            TractDifference {
                locus: 0,
                read: 0,
                offset: 5,
                base: ObservedAllele::A,
            },
        ];
        assert_ne!(shared, scattered);
        let reads_carrying = |list: &[TractDifference]| {
            let mut reads: Vec<u8> = list.iter().map(|d| d.read).collect();
            reads.sort_unstable();
            reads.dedup();
            reads.len()
        };
        assert_eq!(reads_carrying(&shared), 2);
        assert_eq!(reads_carrying(&scattered), 1);
    }

    #[test]
    fn a_flank_difference_and_an_interior_one_come_back_apart() {
        let tract_length = 12_i16;
        let flank = TractDifference {
            locus: 0,
            read: 0,
            offset: -2,
            base: ObservedAllele::T,
        };
        let interior = TractDifference {
            locus: 0,
            read: 0,
            offset: 4,
            base: ObservedAllele::T,
        };
        assert!(flank.offset < 0);
        assert!((0..tract_length).contains(&interior.offset));
    }

    #[test]
    fn the_guard_threshold_is_one_in_ten_of_the_reads_that_differ() {
        let mut records = SsrEvidence::never_walked(2);
        records.walked[0] = true;
        records.offsets[0].add(1, 18);
        records.guard.push(GuardObservation {
            locus: 0,
            length_difference: 1,
            reads: 1,
        });
        assert!(!records.guard_is_over_threshold(0), "1 in 19 is under");

        records.guard.push(GuardObservation {
            locus: 0,
            length_difference: 1,
            reads: 2,
        });
        assert!(records.guard_is_over_threshold(0), "3 in 21 is over");
    }

    // ---- the thirteen values -------------------------------------------------

    fn terms() -> RecordingTerms {
        RecordingTerms {
            selection: selection_terms(),
            kept_loci: CensusLociDigester::new().finish(),
            ssr_stratum_counts: StratumCounts::default(),
            read_cap: ReadCap(100),
            depth_ladder: DepthLadderDigest::of(&DepthBinEdges::new()),
            depth_cap: DepthCap(124),
            coverage_window: Some(Bp(500)),
        }
    }

    fn selection_terms() -> SelectionTerms {
        use crate::ng::parameter_estimation::joint::loci::{
            CatalogBuildSettings, ReferenceDigest, RegionSetDigest,
        };
        use crate::ng::repeat_catalog::StrRepeatCriteria;
        use crate::ng::tandem_repeat::ScanParams;
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
    fn each_of_the_five_unit_values_refuses_on_its_own() {
        // These are the five that an earlier version of the check did not have: two samples
        // can agree on every locus, every seed and every digest and still have written down
        // rows that mean different things.
        let base = terms();

        let mut ladder = base.clone();
        ladder.depth_ladder = DepthLadderDigest([0; 16]);
        assert_eq!(base.first_disagreement(&ladder), Some("depth ladder edges"));

        let mut depth_cap = base.clone();
        depth_cap.depth_cap = DepthCap(60);
        assert_eq!(
            base.first_disagreement(&depth_cap),
            Some("per-position depth cap")
        );

        let mut read_cap = base.clone();
        read_cap.read_cap = ReadCap(30);
        assert_eq!(
            base.first_disagreement(&read_cap),
            Some("per-locus read cap")
        );

        let mut window = base.clone();
        window.coverage_window = Some(Bp(1_000));
        assert_eq!(
            base.first_disagreement(&window),
            Some("coverage window size")
        );

        let mut strata = base.clone();
        let mut counts = StratumCounts::default();
        counts.count(2, 10);
        strata.ssr_stratum_counts = counts;
        assert_eq!(
            base.first_disagreement(&strata),
            Some("per-stratum locus counts")
        );

        assert_eq!(base.first_disagreement(&base.clone()), None);
    }

    #[test]
    fn a_selection_difference_is_named_by_the_selection_rather_than_swallowed() {
        let base = terms();
        let mut other = base.clone();
        other.selection.seed += 1;
        assert_eq!(base.first_disagreement(&other), Some("selection seed"));
    }

    /// **The one value that says what a run *produced* rather than what it was asked for**,
    /// and the only one of them with nothing behind it until now: every fixture in the
    /// module built its terms with the same empty digest, so disabling this comparison
    /// altogether left the whole joint module's suite green. Two samples pooled across it
    /// are two samples whose rows are indexed against different lists of positions.
    #[test]
    fn a_different_set_of_kept_loci_is_refused_and_named() {
        let base = terms();
        let mut other = base.clone();
        let mut digester = CensusLociDigester::new();
        digester.observe(
            0,
            GenomePosition {
                contig: ContigId(0),
                position: Position(1),
            },
        );
        other.kept_loci = digester.finish();
        assert_ne!(base.kept_loci, other.kept_loci);
        assert_eq!(
            base.first_disagreement(&other),
            Some("the loci actually kept")
        );
    }

    // ---- the per-position depth cap ------------------------------------------

    /// **The one function whose stated contract is byte-identity**, and it had no test: a
    /// region-sharded walk and a single-threaded one must thin a deep position's counts the
    /// same way, which is why the share is computed rather than drawn. Turning the rounding
    /// from nearest to downwards leaves every other test in the module passing.
    #[test]
    fn a_thinned_share_rounds_to_nearest_and_never_loses_the_last_read() {
        assert_eq!(thin_to_cap(7, 20, 20), 7, "at the cap, nothing is thinned");
        assert_eq!(
            thin_to_cap(7, 10, 20),
            7,
            "below the cap, nothing is thinned"
        );
        assert_eq!(thin_to_cap(0, 0, 20), 0, "no depth, no reads");
        assert_eq!(
            thin_to_cap(20, 40, 20),
            10,
            "half the depth keeps half the reads, so the fraction survives"
        );
        assert_eq!(thin_to_cap(3, 8, 5), 2, "1.875 rounds to 2, not down to 1");
        assert_eq!(thin_to_cap(2, 9, 5), 1, "1.111 rounds to 1, not up to 2");
        assert_eq!(
            thin_to_cap(1, 400, 20),
            1,
            "a single stray read at 400x is kept rather than rounded away — those reads \
             are what the error rate is fitted from"
        );
        assert_eq!(
            thin_to_cap(400, 400, 20),
            20,
            "and no thinned count ever exceeds the cap"
        );
    }

    // ---- the coverage summary ------------------------------------------------

    fn coverage() -> CoverageByWindow {
        // Ten windows at about three reads a position, one of them short.
        let means = [3.0_f32, 3.0, 3.0, 3.0, 6.0, 6.0, 3.0, 3.0, 3.0, 3.0];
        let positions = [500_u16, 500, 500, 500, 500, 500, 500, 500, 500, 100];
        CoverageByWindow::new(Bp(500), 3.0, &means, &positions, vec![1.0; 20])
    }

    #[test]
    fn a_window_mean_survives_the_byte_it_is_stored_in() {
        let summary = coverage();
        assert!((summary.mean_depth(0) - 3.0).abs() < 0.05);
        // The value the class is about: twice the sample's own median, at three reads a
        // position, where an integer byte could not tell 3 from 4.
        assert!((summary.mean_depth(4) - 6.0).abs() < 0.05);
    }

    #[test]
    fn summing_windows_weights_each_by_its_own_position_count() {
        let summary = coverage();
        // Windows 8 and 9: 500 positions at 3.0 and 100 at 3.0 — the same mean either way.
        assert!((summary.mean_depth_over(8..10) - 3.0).abs() < 0.05);

        // Windows 3..6: 500 at 3.0, 500 at 6.0, 500 at 6.0 → 5.0.
        assert!((summary.mean_depth_over(3..6) - 5.0).abs() < 0.05);

        // The short window matters: a naive mean of windows 0 and 9 weights the 100-position
        // one as if it were full. Here it is not, and the count is what says so.
        assert_eq!(summary.positions(9), 100);
        assert_eq!(summary.positions(0), 500);
    }

    // ---- the writer ----------------------------------------------------------

    use crate::ng::locus_generation::SequenceObservation;

    fn observation(bases: &[u8], group: u32, reads: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.to_vec().into_boxed_slice(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(group),
            num_obs: reads,
            num_fwd: reads,
            q_sum: 0.0,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn generic_locus(
        position: u64,
        reference: &[u8],
        observations: Vec<SequenceObservation>,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(position),
                end: Position(position + reference.len() as u64 - 1),
            },
            reference_bases: reference.to_vec().into_boxed_slice(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Generic,
        }
    }

    /// A selection of three positions on one contig, and a writer over it.
    fn writer_over(positions: &[u64], groups: &[u32]) -> CensusWriter {
        let loci = CensusLoci::from_parts(
            positions
                .iter()
                .map(|p| GenomePosition {
                    contig: ContigId(0),
                    position: Position(*p),
                })
                .collect(),
            Default::default(),
            StratumCounts::default(),
        );
        CensusWriter::new(
            "sample".to_string(),
            &loci,
            groups.iter().map(|g| ReadGroupId(*g)).collect(),
            &|_| Some(ContigId(0)),
            selection_terms(),
            DepthBinEdges::new(),
            ReadCap(100),
            DepthCap(124),
        )
    }

    #[test]
    fn a_walked_position_no_read_reached_is_zero_depth_and_not_never_walked() {
        // Three kept positions; the walk reaches the first two and stops.
        let mut writer = writer_over(&[10, 20, 30], &[0]);
        writer.add_locus(&generic_locus(10, b"A", vec![observation(b"A", 0, 4)]));
        writer.add_locus(&generic_locus(20, b"C", Vec::new()));
        let records = writer.finish(None);

        let generic = &records.generic[&ReadGroupId(0)];
        let edges = DepthBinEdges::new();
        assert_eq!(generic.at(0).0, DepthCode::Binned(edges.bin_for(4)));
        assert_eq!(
            generic.at(1).0,
            DepthCode::Binned(edges.bin_for(0)),
            "walked and empty is data"
        );
        assert_eq!(
            generic.at(2).0,
            DepthCode::NeverWalked,
            "never reached is a bug, and must not read as zero depth"
        );
    }

    /// **The generic locus generator emits nothing at a position no read reached**, so on real
    /// data the previous test's middle case never arises: an uncovered position keeps the
    /// never-walked sentinel, which is the code for a bug. Measured on tomato SRR7279482 at
    /// 25× that is 93,150 kept positions in 1,999,404. Marking the walked stretch is what
    /// separates the two again.
    #[test]
    fn marking_a_walked_stretch_turns_silence_into_zero_depth() {
        let mut writer = writer_over(&[10, 20, 30], &[0]);
        writer.mark_walked(GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(25),
        });
        writer.add_locus(&generic_locus(10, b"A", vec![observation(b"A", 0, 4)]));
        let records = writer.finish(None);

        let generic = &records.generic[&ReadGroupId(0)];
        let edges = DepthBinEdges::new();
        assert_eq!(
            generic.at(0).0,
            DepthCode::Binned(edges.bin_for(4)),
            "a real depth is never overwritten by the mark"
        );
        assert_eq!(
            generic.at(1).0,
            DepthCode::Binned(edges.bin_for(0)),
            "inside the walked stretch and no read reached it: data"
        );
        assert_eq!(
            generic.at(2).0,
            DepthCode::NeverWalked,
            "outside every walked stretch: still a bug"
        );
    }

    #[test]
    fn marking_after_the_loci_arrive_gives_the_same_answer() {
        let mut writer = writer_over(&[10, 20, 30], &[0]);
        writer.add_locus(&generic_locus(10, b"A", vec![observation(b"A", 0, 4)]));
        writer.mark_walked(GenomeRegion {
            contig: ContigId(0),
            start: Position(1),
            end: Position(25),
        });
        let records = writer.finish(None);
        let generic = &records.generic[&ReadGroupId(0)];
        assert_eq!(
            generic.at(0).0,
            DepthCode::Binned(DepthBinEdges::new().bin_for(4))
        );
        assert_eq!(
            generic.at(1).0,
            DepthCode::Binned(DepthBinEdges::new().bin_for(0))
        );
    }

    // ---- the writer at a repeat tract ---------------------------------------

    /// A twelve-base `AT` tract, kept, with the writer over it.
    fn ssr_writer() -> CensusWriter {
        use crate::ng::region_typing::segment_criteria::SsrSegment;
        use crate::ng::repeat_catalog::strata::StratumSampler;

        let segment = SsrSegment::new(
            "c0".into(),
            100,
            111,
            crate::ng::types::Motif::new(b"AT").expect("a two-base motif"),
            1.0,
        )
        .expect("a well-formed tract");
        let mut sampler = StratumSampler::new(16, 0);
        sampler.offer(2, 6, segment);
        let (counts, sample) = sampler.finish();
        let loci = CensusLoci::from_parts(Vec::new(), sample, counts);
        CensusWriter::new(
            "sample".to_string(),
            &loci,
            vec![ReadGroupId(0)],
            &|_| Some(ContigId(0)),
            selection_terms(),
            DepthBinEdges::new(),
            ReadCap(100),
            DepthCap(124),
        )
    }

    fn ssr_locus(observations: Vec<SequenceObservation>) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(100),
                end: Position(111),
            },
            reference_bases: b"ATATATATATAT".to_vec().into_boxed_slice(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(crate::ng::locus_generation::SsrDetail {
                motif: crate::ng::types::Motif::new(b"AT").expect("a two-base motif"),
                left_flank: Box::new([]),
                right_flank: Box::new([]),
            }),
        }
    }

    /// Spec §7.3: two reads carrying the **same** interior substitution must come back as two
    /// entries at one offset — not one entry, and not a count of two. The walk has already
    /// folded identical reads into one observation with a count, so the writer has to unfold
    /// them again or the channel says nothing it was built to say.
    #[test]
    fn two_reads_carrying_one_interruption_are_two_entries_at_one_offset() {
        let mut writer = ssr_writer();
        writer.add_locus(&ssr_locus(vec![observation(b"ATATGTATATAT", 0, 2)]));
        let records = writer.finish(None);
        let ssr = &records.ssr[&ReadGroupId(0)];

        assert_eq!(
            ssr.differences().len(),
            2,
            "two reads, one interruption each"
        );
        assert_eq!(ssr.differences()[0].offset, 4);
        assert_eq!(ssr.differences()[1].offset, 4);
        assert_ne!(
            ssr.differences()[0].read,
            ssr.differences()[1].read,
            "an allele on two reads, not one read seen twice"
        );
        assert_eq!(ssr.differences()[0].base, ObservedAllele::G);
    }

    /// A read whose tract slipped a whole unit has no base-for-base correspondence with the
    /// reference, so it contributes its offset and **nothing to the denominator**. If it did,
    /// the STR error rate would be a ratio of two quantities counted over different reads.
    #[test]
    fn a_slipped_read_contributes_a_length_and_no_base_comparison() {
        let mut writer = ssr_writer();
        writer.add_locus(&ssr_locus(vec![
            observation(b"ATATATATAT", 0, 3),
            observation(b"ATATATATATAT", 0, 5),
        ]));
        let records = writer.finish(None);
        let ssr = &records.ssr[&ReadGroupId(0)];

        assert_eq!(ssr.offsets(0).at(-1), 3, "one unit short");
        assert_eq!(ssr.offsets(0).at(0), 5);
        assert_eq!(
            ssr.bases_compared(0),
            5 * 12,
            "only the five reads at the reference length were compared"
        );
        assert!(ssr.differences().is_empty());
    }

    #[test]
    fn a_read_group_that_saw_nothing_here_still_gets_its_zero() {
        // Group 1 puts no read at this position. Its record must say zero, not never walked
        // — the entry is the denominator its own error rate is fitted against.
        let mut writer = writer_over(&[10], &[0, 1]);
        writer.add_locus(&generic_locus(10, b"A", vec![observation(b"A", 0, 6)]));
        let records = writer.finish(None);

        let edges = DepthBinEdges::new();
        assert_eq!(
            records.generic[&ReadGroupId(0)].at(0).0,
            DepthCode::Binned(edges.bin_for(6))
        );
        assert_eq!(
            records.generic[&ReadGroupId(1)].at(0).0,
            DepthCode::Binned(edges.bin_for(0))
        );
    }

    #[test]
    fn depth_is_the_read_groups_own_and_not_the_samples() {
        // Pooling would give both groups a depth of 7; each saw its own.
        let mut writer = writer_over(&[10], &[0, 1]);
        writer.add_locus(&generic_locus(
            10,
            b"A",
            vec![observation(b"A", 0, 5), observation(b"A", 1, 2)],
        ));
        let records = writer.finish(None);

        let edges = DepthBinEdges::new();
        assert_eq!(
            records.generic[&ReadGroupId(0)].at(0).0,
            DepthCode::Binned(edges.bin_for(5))
        );
        assert_eq!(
            records.generic[&ReadGroupId(1)].at(0).0,
            DepthCode::Binned(edges.bin_for(2))
        );
    }

    #[test]
    fn a_non_reference_read_lands_on_its_own_base_and_a_matching_one_on_none() {
        let mut writer = writer_over(&[10, 20], &[0]);
        writer.add_locus(&generic_locus(
            10,
            b"A",
            vec![observation(b"A", 0, 3), observation(b"G", 0, 2)],
        ));
        writer.add_locus(&generic_locus(20, b"C", vec![observation(b"C", 0, 4)]));
        let records = writer.finish(None);

        let generic = &records.generic[&ReadGroupId(0)];
        let (depth, alleles) = generic.at(0);
        assert_eq!(depth, DepthCode::Binned(DepthBinEdges::new().bin_for(5)));
        assert_eq!(
            alleles,
            [AlleleObservation {
                index: 0,
                allele: ObservedAllele::G,
                reads: 2,
            }]
        );
        assert!(
            generic.at(1).1.is_empty(),
            "a position whose reads all matched carries no sparse entry"
        );
    }

    #[test]
    fn a_position_the_selection_does_not_hold_is_ignored() {
        let mut writer = writer_over(&[10], &[0]);
        writer.add_locus(&generic_locus(11, b"A", vec![observation(b"T", 0, 9)]));
        let records = writer.finish(None);
        assert_eq!(
            records.generic[&ReadGroupId(0)].at(0).0,
            DepthCode::NeverWalked
        );
    }

    #[test]
    fn the_kept_loci_digest_covers_every_kept_locus_even_where_the_walk_stopped() {
        // Two writers over the same selection, one of which walked less. The digest witnesses
        // the selection, so they agree — and a third over a different selection does not.
        let mut walked = writer_over(&[10, 20, 30], &[0]);
        walked.add_locus(&generic_locus(10, b"A", vec![observation(b"A", 0, 3)]));
        let mut idle = writer_over(&[10, 20, 30], &[0]);
        idle.add_locus(&generic_locus(30, b"A", vec![observation(b"A", 0, 3)]));
        assert_eq!(
            walked.finish(None).terms.kept_loci,
            idle.finish(None).terms.kept_loci
        );

        let other = writer_over(&[10, 20, 31], &[0]);
        assert_ne!(
            writer_over(&[10, 20, 30], &[0])
                .finish(None)
                .terms
                .kept_loci,
            other.finish(None).terms.kept_loci
        );
    }

    #[test]
    fn how_many_windows_to_sum_falls_out_of_the_samples_own_depth() {
        // 25× at 500 bp is 12,500 aligned bases — one window is enough.
        let deep = CoverageByWindow::new(Bp(500), 25.0, &[25.0], &[500], Vec::new());
        assert_eq!(deep.windows_to_sum(), 1);

        // 2.5× at 500 bp is 1,250 — ten windows, which is the 5 kb the measurement found.
        let shallow = CoverageByWindow::new(Bp(500), 2.5, &[2.5], &[500], Vec::new());
        assert_eq!(shallow.windows_to_sum(), 10);
        assert_eq!(
            shallow.windows_to_sum() * 500,
            5_000,
            "5 kb is what 2.5x needs, measured on eight tomato samples"
        );
    }
}
