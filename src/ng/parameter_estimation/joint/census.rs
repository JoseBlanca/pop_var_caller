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
use std::ops::RangeInclusive;

use md5::{Digest, Md5};

use crate::ng::locus_generation::{LocusKind, ReadWitness, SampleLocusObservations};
use crate::ng::parameter_estimation::generic::depth_bins::{DepthBin, DepthBinEdges};
use crate::ng::parameter_estimation::joint::loci::{
    CensusLoci, CensusLociDigest, CensusLociDigester, SelectionTerms,
};
use crate::ng::repeat_catalog::StratumCounts;
use crate::ng::types::{ContigId, GenomePosition, GenomeRegion, Position, ReadGroupId};

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

/// How many codes fit in one entry. Thirty bins plus the sentinel is 31, so five bits.
pub const DEPTH_CODE_BITS: u32 = 5;

/// The value [`DepthCode::NeverWalked`] is stored as — the top of the five-bit range, so
/// adding a rung to the ladder collides with it loudly rather than shifting it.
const NEVER_WALKED_CODE: u8 = (1 << DEPTH_CODE_BITS) - 1;

// **The ladder cannot outgrow the field it is stored in without failing to compile.** The
// bins take codes `0..bin_count` and the sentinel takes the top of the range, so the last
// bin has to sit strictly below it. `to_bits` asserts the same thing per entry, and that
// assert fires while a run is walking a genome; this one fires while it is being built.
const _: () = assert!(
    crate::ng::parameter_estimation::generic::depth_bins::CENSUS_DEPTH_BIN_COUNT
        <= NEVER_WALKED_CODE as usize,
    "the census ladder has outgrown the five-bit depth code: its top rung would be written \
     as the never-walked sentinel, which is the code for a bug"
);

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
    /// **Exact, not binned, and one byte.** Nearly every entry is a miscall of one to three
    /// reads, so the tail a ladder would compress is empty at any depth and binning would
    /// buy a fourteenth bin width in [`RecordingTerms`] for nothing.
    ///
    /// **One byte is safe because [`DepthCap`] bounds it exactly**: a count cannot exceed
    /// its position's depth, the walk thins a position's counts by the same ratio it thins
    /// the depth, and the cap refuses a value a byte cannot hold. An entry is an index, an
    /// allele and a count — six bytes here against about nine with a wider count, which at
    /// 100× is a megabyte a read group.
    pub reads: u8,
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
    /// How far into the tract the mismatching base sat, in bases from its first: `0..len`,
    /// and nothing else.
    ///
    /// **The sequence either side of the tract is not recorded here — DECIDED 2026-08-14
    /// (owner).** It is present in the locus only because the aligner needs it to anchor a
    /// read, and it is ordinary non-repetitive sequence: what happens there is the generic
    /// half of the census's business, at the positions it already keeps. So the STR half
    /// compares a read against the tract and against nothing else, and an earlier version of
    /// this comment promising negative offsets described a state no code path produced.
    ///
    /// Signed because reading it as a plain index would hide a later change of mind.
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
///
/// **It thins the allele counts and no longer clips the depth — CHANGED 2026-08-14.** The
/// depth code stores what the position actually held; the counts beside it were thinned to
/// this cap, keeping the fractions they showed. So a consumer asking *how many copies is
/// this* reads the depth, and one dividing a count by a depth must first put the depth into
/// the counts' own units — which is [`denominator_for`](Self::denominator_for), and nothing
/// else records it.
///
/// **It refuses a value a byte cannot hold** ([`new`](Self::new)), because it is the bound on
/// [`AlleleObservation::reads`]: the counts are thinned to it, so a cap above 255 would let
/// them saturate silently while the depth field beside them said otherwise. Nothing else
/// stops the two drifting apart — the cap is a run-time value — so the constructor is where
/// the relationship is a property rather than a comment.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct DepthCap(u32);

impl DepthCap {
    /// The largest cap the record's one-byte count can hold, and the value a fit with no
    /// sample at all falls back to.
    pub const MAX: Self = Self(u8::MAX as u32);

    /// A cap of `cap` reads a position.
    ///
    /// # Panics
    ///
    /// Above `u8::MAX`, which is the whole point: [`AlleleObservation::reads`] is one byte
    /// and this cap is what bounds it. `const fn`, so a cap written as a constant is refused
    /// while the run is being built rather than while it is walking a genome — the same move
    /// the ladder's own five-bit assertion makes.
    #[must_use]
    pub const fn new(cap: u32) -> Self {
        assert!(
            cap <= u8::MAX as u32,
            "a per-position depth cap above 255 cannot bound a one-byte allele count: the \
             counts would saturate silently while the depth beside them said otherwise"
        );
        Self(cap)
    }

    /// The depth above which a position's reads are subsampled.
    #[must_use]
    pub const fn get(self) -> u32 {
        self.0
    }

    /// The depths the allele counts at a position were counted out of, given the depths its
    /// stored code stands for.
    ///
    /// **A range in, a range out**, because the stored code is a bin: the count of reads that
    /// disagreed with the reference is exact and the depth is a span, so anything dividing one
    /// by the other has to carry the span or it inherits the bin's width — an error of up to a
    /// sixth at thirty reads a position
    /// (`doc/devel/ng/reports/contamination_floor_and_duplicated_class_2026-08-13.md` §4).
    ///
    /// **Both ends are clamped, which is what collapses a deep position to a point.** At a cap
    /// of 124 a position holding 300 reads has a code standing for 263 to 336; its counts were
    /// thinned to 124, so the denominator is 124 to 124 and the fit sees the one depth the
    /// counts were actually taken out of.
    #[must_use]
    pub fn denominator_for(self, depths: RangeInclusive<u32>) -> RangeInclusive<u32> {
        (*depths.start()).min(self.0)..=(*depths.end()).min(self.0)
    }
}

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

/// The twelve values the fit refuses to pool across.
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
        None
    }
}

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
                counts[observation.allele.code() as usize] += u32::from(observation.reads);
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
///
/// **It returns the byte the record stores**, which is safe for the reason [`DepthCap`]
/// exists: a count cannot exceed its position's depth, the surviving share cannot exceed the
/// cap, and the cap cannot exceed 255. The conversion says so rather than truncating — a
/// truncated count is a wrong allele fraction and has no symptom.
///
/// # Panics
///
/// When the surviving count will not fit in a byte, which means an allele carried more reads
/// than the position it sits at held.
fn thin_to_cap(reads: u32, depth: u32, cap: DepthCap) -> u8 {
    let cap = cap.get();
    let surviving = if depth <= cap || depth == 0 {
        reads
    } else {
        let scaled = f64::from(reads) * f64::from(cap) / f64::from(depth);
        (scaled.round() as u32).min(cap).max(u32::from(reads > 0))
    };
    u8::try_from(surviving).unwrap_or_else(|_| {
        panic!(
            "{surviving} reads on one allele where the position held {depth} and the cap is \
             {cap}: a count cannot exceed the depth it was thinned to, and the record holds \
             one byte"
        )
    })
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
            // **The depth recorded is the one the position actually held** (§2.2, changed
            // 2026-08-14). The cap below still thins the allele counts, which is what keeps
            // them inside their own field and what bounds the sparse list at high depth; it
            // no longer clips this. Nothing is lost by that and one thing is gained: a
            // consumer asking *how many copies of this stretch does the sample carry* reads
            // a depth that can answer, where a clipped one went blind at twice the cap.
            // A consumer needing the counts' own denominator computes it with
            // [`DepthCap::denominator_for`], the cap travelling in the recording terms.
            for (group, depth) in &depths {
                let records = self
                    .generic
                    .entry(*group)
                    .or_insert_with(|| GenericEvidence::never_walked(self.generic_loci.len()));
                records
                    .depth
                    .set(index, DepthCode::Binned(self.edges.bin_for(depth[offset])));
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
                let reads = thin_to_cap(observation.num_obs, group_depth, self.depth_cap);
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

    /// The finished evidence, with the digest of what was actually written.
    pub fn finish(mut self) -> SampleCensusEvidence {
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
        SampleCensusEvidence {
            sample: self.sample,
            generic: self.generic,
            ssr: self.ssr,
            terms: RecordingTerms {
                selection: self.terms,
                kept_loci: self.digester.finish(),
                ssr_stratum_counts: self.stratum_counts,
                read_cap: self.read_cap,
                depth_ladder: DepthLadderDigest::of(&self.edges),
                depth_cap: self.depth_cap,
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
        let edges = DepthBinEdges::for_census();
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
            DepthCode::Binned(DepthBinEdges::for_census().bin_for(0))
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
        let edges = DepthBinEdges::for_census();
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
        let edges = DepthBinEdges::for_census();

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

    /// **The sequence either side of a tract is not this half's business** (owner,
    /// 2026-08-14): it sits in the locus only so the aligner can anchor a read, and the
    /// generic half of the census already keeps the ordinary positions it covers. So every
    /// mismatch recorded here lands inside the tract, at its own distance from the tract's
    /// first base.
    ///
    /// The test this replaced asserted that a hand-written `-2` was below zero and called no
    /// code at all, standing in for a property the writer could not produce.
    #[test]
    fn a_mismatch_is_recorded_at_its_own_distance_into_the_tract_and_never_outside_it() {
        // A twelve-base `AT` tract, and one read whose ninth base is a G.
        let mut writer = ssr_writer();
        writer.add_locus(&ssr_locus(vec![observation(b"ATATATATGTAT", 0, 1)]));
        let evidence = writer.finish();
        let ssr = &evidence.ssr[&ReadGroupId(0)];

        assert_eq!(ssr.differences().len(), 1);
        assert_eq!(
            ssr.differences()[0].offset,
            8,
            "the ninth base of the tract, counted from its first"
        );
        assert_eq!(ssr.differences()[0].base, ObservedAllele::G);
        assert!(
            ssr.differences()
                .iter()
                .all(|difference| (0..12).contains(&difference.offset)),
            "nothing outside the tract is recorded, so no offset is negative or past its end"
        );
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

    // ---- the twelve values ----------------------------------------------

    fn terms() -> RecordingTerms {
        RecordingTerms {
            selection: selection_terms(),
            kept_loci: CensusLociDigester::new().finish(),
            ssr_stratum_counts: StratumCounts::default(),
            read_cap: ReadCap(100),
            depth_ladder: DepthLadderDigest::of(&DepthBinEdges::for_census()),
            depth_cap: DepthCap::new(124),
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
    fn each_of_the_four_unit_values_refuses_on_its_own() {
        // These are the four that an earlier version of the check did not have: two samples
        // can agree on every locus, every seed and every digest and still have written down
        // rows that mean different things.
        let base = terms();

        let mut ladder = base.clone();
        ladder.depth_ladder = DepthLadderDigest([0; 16]);
        assert_eq!(base.first_disagreement(&ladder), Some("depth ladder edges"));

        let mut depth_cap = base.clone();
        depth_cap.depth_cap = DepthCap::new(60);
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

    /// **Spec §7.11.** A position above the cap keeps the depth it actually held, and only its
    /// allele counts are thinned — so the fraction the reads showed survives, and a consumer
    /// recovers the counts' own denominator with `min(depth, cap)`.
    ///
    /// **This is the step's silent failure made loud.** Before 2026-08-14 the depth was
    /// clipped here too, and a consumer that kept reading it that way after the change would
    /// divide a thinned count by an unthinned depth and report a rate several times too low,
    /// with nothing to notice.
    #[test]
    fn a_position_above_the_cap_keeps_its_own_depth_and_thins_only_the_counts() {
        // Forty reads at a cap of twenty: a quarter of them on a non-reference base.
        let mut writer = writer_over_capped(&[10], &[0], DepthCap::new(20));
        writer.add_locus(&generic_locus(
            10,
            b"A",
            vec![observation(b"A", 0, 30), observation(b"G", 0, 10)],
        ));
        let evidence = writer.finish();
        let (depth, alleles) = evidence.generic[&ReadGroupId(0)].at(0);

        let edges = DepthBinEdges::for_census();
        let DepthCode::Binned(bin) = depth else {
            panic!("a walked position is not never-walked")
        };
        assert!(
            edges.depth_range(bin).contains(&40),
            "the code stands for the depth the position actually held, not for the cap"
        );
        assert!(
            !edges.depth_range(bin).contains(&20),
            "and it is not the cap's own bin"
        );

        assert_eq!(alleles.len(), 1);
        assert_eq!(
            alleles[0].reads, 5,
            "a quarter of the cap, as a quarter of the depth was"
        );
        assert!(
            alleles[0].reads <= 20,
            "no count survives above the cap it was thinned to"
        );

        // What a consumer must do, and the only thing that recovers the fraction.
        let denominator = DepthCap::new(20).denominator_for(edges.depth_range(bin));
        assert_eq!(denominator, 20..=20, "min(depth, cap) is the cap here");
        assert!(
            (f64::from(alleles[0].reads) / f64::from(*denominator.start()) - 0.25).abs() < 1e-12,
            "the fraction the reads showed survives the thinning"
        );
    }

    /// Below the cap the two are the same number, which is the case every other fixture in
    /// this module sits in — so this is what says the clamp is not simply always firing.
    #[test]
    fn a_position_below_the_cap_is_its_own_denominator() {
        let mut writer = writer_over_capped(&[10], &[0], DepthCap::new(20));
        writer.add_locus(&generic_locus(
            10,
            b"A",
            vec![observation(b"A", 0, 4), observation(b"G", 0, 2)],
        ));
        let evidence = writer.finish();
        let (depth, alleles) = evidence.generic[&ReadGroupId(0)].at(0);

        let edges = DepthBinEdges::for_census();
        assert_eq!(depth, DepthCode::Binned(edges.bin_for(6)));
        assert_eq!(alleles[0].reads, 2, "nothing is thinned below the cap");
        let DepthCode::Binned(bin) = depth else {
            panic!("a walked position is not never-walked")
        };
        assert_eq!(
            DepthCap::new(20).denominator_for(edges.depth_range(bin)),
            6..=6,
            "below the cap the denominator is the depth itself, exactly"
        );
    }

    /// **The one function whose stated contract is byte-identity**, and it had no test: a
    /// region-sharded walk and a single-threaded one must thin a deep position's counts the
    /// same way, which is why the share is computed rather than drawn. Turning the rounding
    /// from nearest to downwards leaves every other test in the module passing.
    #[test]
    fn a_thinned_share_rounds_to_nearest_and_never_loses_the_last_read() {
        assert_eq!(
            thin_to_cap(7, 20, DepthCap::new(20)),
            7,
            "at the cap, nothing is thinned"
        );
        assert_eq!(
            thin_to_cap(7, 10, DepthCap::new(20)),
            7,
            "below the cap, nothing is thinned"
        );
        assert_eq!(
            thin_to_cap(0, 0, DepthCap::new(20)),
            0,
            "no depth, no reads"
        );
        assert_eq!(
            thin_to_cap(20, 40, DepthCap::new(20)),
            10,
            "half the depth keeps half the reads, so the fraction survives"
        );
        assert_eq!(
            thin_to_cap(3, 8, DepthCap::new(5)),
            2,
            "1.875 rounds to 2, not down to 1"
        );
        assert_eq!(
            thin_to_cap(2, 9, DepthCap::new(5)),
            1,
            "1.111 rounds to 1, not up to 2"
        );
        assert_eq!(
            thin_to_cap(1, 400, DepthCap::new(20)),
            1,
            "a single stray read at 400x is kept rather than rounded away — those reads \
             are what the error rate is fitted from"
        );
        assert_eq!(
            thin_to_cap(400, 400, DepthCap::new(20)),
            20,
            "and no thinned count ever exceeds the cap"
        );
    }

    /// **The refusal that makes the one-byte count safe** (spec §2.1, §2.2). The cap is a
    /// run-time value and the count is a byte, and nothing but this stops the two drifting
    /// apart: a cap of 300 would let a position's counts be thinned to 300 and written as
    /// 44, with the depth beside them saying 300 and no symptom anywhere.
    #[test]
    #[should_panic(expected = "cannot bound a one-byte allele count")]
    fn a_cap_a_byte_cannot_hold_is_refused_at_construction() {
        let _ = DepthCap::new(u8::MAX as u32 + 1);
    }

    #[test]
    fn the_largest_cap_a_byte_holds_is_accepted() {
        // The boundary from the other side: 255 is a cap, 256 is not, and `MAX` is the first.
        assert_eq!(DepthCap::new(255).get(), 255);
        assert_eq!(DepthCap::MAX, DepthCap::new(u8::MAX as u32));
    }

    /// A position sitting exactly at the cap keeps every read it showed: the count that
    /// fills the byte is the one the cap allows, so the widest count the encoding can meet
    /// is written and read back unchanged rather than wrapping to nothing.
    #[test]
    fn a_position_at_the_cap_round_trips_at_the_widest_count_the_byte_holds() {
        let cap = DepthCap::MAX;
        let mut writer = writer_over_capped(&[10], &[0], cap);
        writer.add_locus(&generic_locus(
            10,
            b"A",
            vec![observation(b"G", 0, cap.get())],
        ));
        let evidence = writer.finish();
        let (depth, alleles) = evidence.generic[&ReadGroupId(0)].at(0);

        let edges = DepthBinEdges::for_census();
        let DepthCode::Binned(bin) = depth else {
            panic!("a walked position is not never-walked")
        };
        assert!(
            edges.depth_range(bin).contains(&cap.get()),
            "the position held exactly the cap, and its code says so"
        );
        assert_eq!(
            alleles[0].reads,
            u8::MAX,
            "255 reads on one allele is the widest count a byte holds, and it survives whole"
        );
        assert_eq!(
            evidence.allele_counts(0)[ObservedAllele::G.code() as usize],
            u32::from(u8::MAX),
            "and it folds across read groups without being truncated on the way"
        );
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
        writer_over_capped(positions, groups, DepthCap::new(124))
    }

    /// The same, at a chosen per-position cap. **A cap below the ladder's ceiling is the only
    /// way to see the cap at all**: it and the ladder's top used to be one number, and the
    /// tests that came with them all sat at depths neither reaches.
    fn writer_over_capped(positions: &[u64], groups: &[u32], cap: DepthCap) -> CensusWriter {
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
            DepthBinEdges::for_census(),
            ReadCap(100),
            cap,
        )
    }

    #[test]
    fn a_walked_position_no_read_reached_is_zero_depth_and_not_never_walked() {
        // Three kept positions; the walk reaches the first two and stops.
        let mut writer = writer_over(&[10, 20, 30], &[0]);
        writer.add_locus(&generic_locus(10, b"A", vec![observation(b"A", 0, 4)]));
        writer.add_locus(&generic_locus(20, b"C", Vec::new()));
        let records = writer.finish();

        let generic = &records.generic[&ReadGroupId(0)];
        let edges = DepthBinEdges::for_census();
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
        let records = writer.finish();

        let generic = &records.generic[&ReadGroupId(0)];
        let edges = DepthBinEdges::for_census();
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
        let records = writer.finish();
        let generic = &records.generic[&ReadGroupId(0)];
        assert_eq!(
            generic.at(0).0,
            DepthCode::Binned(DepthBinEdges::for_census().bin_for(4))
        );
        assert_eq!(
            generic.at(1).0,
            DepthCode::Binned(DepthBinEdges::for_census().bin_for(0))
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
            DepthBinEdges::for_census(),
            ReadCap(100),
            DepthCap::new(124),
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
        let records = writer.finish();
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
        let records = writer.finish();
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
        let records = writer.finish();

        let edges = DepthBinEdges::for_census();
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
        let records = writer.finish();

        let edges = DepthBinEdges::for_census();
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
        let records = writer.finish();

        let generic = &records.generic[&ReadGroupId(0)];
        let (depth, alleles) = generic.at(0);
        assert_eq!(
            depth,
            DepthCode::Binned(DepthBinEdges::for_census().bin_for(5))
        );
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
        let records = writer.finish();
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
            walked.finish().terms.kept_loci,
            idle.finish().terms.kept_loci
        );

        let other = writer_over(&[10, 20, 31], &[0]);
        assert_ne!(
            writer_over(&[10, 20, 30], &[0]).finish().terms.kept_loci,
            other.finish().terms.kept_loci
        );
    }
}
