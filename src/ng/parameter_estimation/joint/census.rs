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
//! walked with no coverage (data), and a locus whose reads all matched (data). [`DepthCode`]
//! is what makes the first expressible at an ordinary position, and [`WalkedBits`] at a repeat
//! tract.
//!
//! **A repeat tract had a fourth until 2026-08-14** — reads reached it and crossed none of it,
//! so none reported a length. That count is now kept per *stratum* rather than per locus
//! ([`SsrEvidence::covering_not_crossing`]), which is the grain the loss runs along: a tract
//! longer than the reads is never crossed in any sample at any depth.
//!
//! # What is not here yet
//!
//! **The byte-level census file.** Where the evidence lives between the walk and the fit is a
//! non-goal of the spec (§1.2) — what it requires is a property, that the fit reaches every
//! sample's evidence without walking the reads again. The encoding that carries the *content*
//! is here and is tested by packing and unpacking it; the framing that would put it in a file
//! is the next unit, and it changes none of these types.

use std::collections::BTreeMap;
use std::io::{Read, Seek};
use std::ops::RangeInclusive;
use std::path::PathBuf;

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

    /// The allele a stored code stands for — `None` for a code no allele has, which is what a
    /// census file written by something else looks like.
    pub fn of_code(code: u8) -> Option<Self> {
        match code {
            0 => Some(Self::A),
            1 => Some(Self::C),
            2 => Some(Self::G),
            3 => Some(Self::T),
            4 => Some(Self::Other),
            _ => None,
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

    /// The array a census file wrote down, read back.
    ///
    /// # Panics
    ///
    /// When the bytes are not the number `len` entries at five bits each need. A short array
    /// would read a position off the end of the buffer, and a long one means the writer and the
    /// reader disagree about the encoding — neither has a symptom in the values.
    pub fn from_bytes(bits: Vec<u8>, len: usize) -> Self {
        assert_eq!(
            bits.len(),
            len.saturating_mul(DEPTH_CODE_BITS as usize).div_ceil(8),
            "{len} depth codes at {DEPTH_CODE_BITS} bits each are not {} bytes",
            bits.len()
        );
        Self { bits, len }
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

    /// The buckets themselves, `-4 … +4` in order — for the codec that writes them.
    pub fn counts(&self) -> &[u16; OFFSET_BUCKETS] {
        &self.counts
    }

    /// The buckets a census file wrote down, read back.
    pub fn from_counts(counts: [u16; OFFSET_BUCKETS]) -> Self {
        Self { counts }
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
    /// Which tract, **numbered within its own stratum** rather than across the whole STR
    /// selection: a section holds one stratum's tracts, so a genome-wide index would name a
    /// tract the section does not have.
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
    /// Which tract, **numbered within its own stratum** rather than across the whole STR
    /// selection — a section holds one stratum's tracts (see [`SsrEvidence`]).
    pub locus: u32,
    /// Which of this locus's reads carried it, numbered from zero within the locus and the
    /// read group, over the reads whose tract was the reference's length — the ones a
    /// mismatch can be read off at all.
    ///
    /// **Two entries at one offset on one read is a different observation from the same two
    /// on two reads**, and a read-blind encoding passes every other check in the suite. The
    /// numbering therefore has to run across the whole locus: the walk folds reads carrying
    /// identical bases into one observation, so a counter restarting at each observation gave
    /// two reads that differ in different places the same number.
    ///
    /// **Two bytes, not one — WIDENED 2026-08-14 (owner).** A byte stops at 255 where a locus
    /// enters up to 1,000 reads ([`ReadCap`]), and on the human benchmark trio — about 265
    /// crossing reads a tract a sample — **one mismatch in five came back numbered 255**, each
    /// one a read that could no longer be told from the others there. Two bytes reach 65,535,
    /// which is past any cap this record is written under, so the numbering cannot collapse at
    /// the deep end of the range the caller commits to. It costs four bytes an entry, on a list
    /// the error rate fills rather than the variants.
    pub read: u16,
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

/// Which stratum a tract is in: its motif length and the reference's repeat count.
///
/// **It is what two of a tract's counts are held per, rather than per locus** — see
/// [`SsrEvidence`] — and it is what the slippage numbers are fitted within, so the two modules
/// name one type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Stratum {
    /// How many bases one repeat unit is.
    pub period: u8,
    /// How many copies of it the reference carries.
    pub reference_repeats: u64,
}

/// One bit a locus: whether the walk reached it.
///
/// **A byte would carry 58 kB of information in 0.46 MB** on tomato's 462,701 kept tracts a
/// read group, and this is the field that has to exist at every locus whether or not anything
/// happened there — it is the STR half's answer to the generic half's never-walked sentinel,
/// and every other vector reads the same for a locus never walked and one walked with no read.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct WalkedBits {
    bits: Vec<u8>,
    len: usize,
}

impl WalkedBits {
    /// `len` loci, none of them walked.
    pub fn none_of(len: usize) -> Self {
        Self {
            bits: vec![0; len.div_ceil(8)],
            len,
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// # Panics
    ///
    /// When `locus` is past the end — the evidence is sized from the selection, so an index
    /// outside it means the writer and the selection disagree.
    pub fn set(&mut self, locus: usize) {
        assert!(locus < self.len, "locus {locus} is outside the selection");
        self.bits[locus / 8] |= 1 << (locus % 8);
    }

    pub fn get(&self, locus: usize) -> bool {
        assert!(locus < self.len, "locus {locus} is outside the selection");
        self.bits[locus / 8] & (1 << (locus % 8)) != 0
    }

    /// The packed bytes — one eighth of what a `bool` a locus costs.
    pub fn as_bytes(&self) -> &[u8] {
        &self.bits
    }

    /// The bits a census file wrote down, read back.
    ///
    /// # Panics
    ///
    /// When the bytes are not the number `len` loci need — see
    /// [`PackedDepthCodes::from_bytes`].
    pub fn from_bytes(bits: Vec<u8>, len: usize) -> Self {
        assert_eq!(
            bits.len(),
            len.div_ceil(8),
            "{len} loci at one bit each are not {} bytes",
            bits.len()
        );
        Self { bits, len }
    }
}

/// One read group's tracts **for one stratum** — the smallest piece of the STR half anything
/// asks for, and what [`SectionKey::Ssr`] names.
///
/// **One stratum and not a whole read group — since 2026-08-14.** The slippage numbers are
/// fitted within a stratum and the fit reads one band of strata at a time, so a value spanning
/// every stratum is a value no caller can use a part of
/// (`parameter_prepass_joint_records.md` §6.2). Every per-tract vector here is therefore indexed
/// by the tract's position **within this stratum**, in genome order, and the two counts that
/// spanned strata are now one number each.
///
/// **Three states at a locus**: the region was never walked; it was walked and no read crossed
/// the whole tract, so none reported a length; or reads crossed it. The generic set's
/// zero-depth-against-quiet distinction is the last pair, and [`walked`](WalkedBits) is what
/// makes the first expressible.
///
/// **A fourth state was expressible per locus until 2026-08-14 and is not any more** — reads
/// that reached the tract and crossed none of it. Spec §3 moved that count to the stratum, so
/// what it distinguishes now is a *stratum* whose tracts are longer than the reads from one
/// that was unlucky with coverage, which is the axis the censoring actually runs along and the
/// grain the count is reported at. Per locus the fit never acted on it: it summed the count
/// the moment it read it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SsrEvidence {
    offsets: Vec<OffsetCounts>,
    /// This stratum's censored reads, over every tract in it (spec §3).
    covering_not_crossing: u64,
    walked: WalkedBits,
    /// This stratum's base-comparison denominator (spec §3): the rate is fitted per stratum,
    /// and per locus it was 1.85 MB a read group on tomato that nothing read.
    bases_compared: u64,
    guard: Vec<GuardObservation>,
    differences: Vec<TractDifference>,
}

impl SsrEvidence {
    /// A stratum's `loci` tracts, none of them walked.
    pub fn never_walked(loci: usize) -> Self {
        Self {
            offsets: vec![OffsetCounts::default(); loci],
            covering_not_crossing: 0,
            walked: WalkedBits::none_of(loci),
            bases_compared: 0,
            guard: Vec::new(),
            differences: Vec::new(),
        }
    }

    /// A section a census file wrote down, read back.
    ///
    /// # Panics
    ///
    /// When the offsets and the walked bits are not the same length, or when an entry names a
    /// tract this stratum does not have. Both leave a consumer reading one tract's evidence at
    /// another's, and neither has a symptom.
    #[allow(
        clippy::too_many_arguments,
        reason = "every field of the section, and a struct to bundle them would be the section"
    )]
    pub fn from_parts(
        offsets: Vec<OffsetCounts>,
        covering_not_crossing: u64,
        walked: WalkedBits,
        bases_compared: u64,
        guard: Vec<GuardObservation>,
        differences: Vec<TractDifference>,
    ) -> Self {
        assert_eq!(
            offsets.len(),
            walked.len(),
            "a section's offsets and its walked bits count different tracts"
        );
        let loci = offsets.len() as u32;
        assert!(
            guard.iter().all(|entry| entry.locus < loci)
                && differences.iter().all(|entry| entry.locus < loci),
            "an entry names a tract this stratum does not hold"
        );
        Self {
            offsets,
            covering_not_crossing,
            walked,
            bases_compared,
            guard,
            differences,
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

    /// Reads that reached a tract of this stratum and crossed no whole copy of it — a censored
    /// lower bound on the length, and **counted per stratum**.
    ///
    /// **The censoring is not random**, which is why it is recorded at all: a tract longer than
    /// a read is never crossed, in every sample at every depth, so the loss runs along repeat
    /// count — which is what a stratum is — and a stratum unreadable at this read length must
    /// not arrive looking like one that was merely unlucky with coverage. On the human
    /// benchmark trio 37 reads in every 100 that reach a tract never cross it, and on the
    /// 63-accession tomato cohort it is close to half.
    ///
    /// **Per locus it was 0.93 MB a read group that nothing read**: the fit added it into a
    /// running per-stratum total the moment it saw it (spec §3).
    pub fn covering_not_crossing(&self) -> u64 {
        self.covering_not_crossing
    }

    /// The denominator the STR substitution rate is fitted against — one number, because this
    /// section *is* one stratum and that is the grain the rate is fitted at.
    pub fn bases_compared(&self) -> u64 {
        self.bases_compared
    }

    pub fn guard(&self) -> &[GuardObservation] {
        &self.guard
    }

    pub fn differences(&self) -> &[TractDifference] {
        &self.differences
    }

    /// Which of the three states this locus is in.
    pub fn state(&self, locus: usize) -> SsrLocusState {
        if !self.walked.get(locus) {
            SsrLocusState::NeverWalked
        } else if self.offsets[locus].total() > 0 {
            SsrLocusState::Crossed
        } else {
            SsrLocusState::NoLengthReported
        }
    }

    /// Whether the walk reached this locus at all.
    pub fn walked(&self, locus: usize) -> bool {
        self.walked.get(locus)
    }

    /// The walked bits themselves, for a consumer measuring what the evidence costs.
    pub fn walked_bits(&self) -> &WalkedBits {
        &self.walked
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

/// Which of the three states a sample is in at a kept STR locus.
///
/// **There were four until 2026-08-14**, the fourth being *reads reached the locus and crossed
/// none of it*. That count moved to the stratum (spec §3), so per locus it can no longer be
/// told from *no read reached the locus at all* — both report no length, which is the only
/// thing the fit does with either. The distinction survives where it is acted on: per stratum,
/// in [`SsrEvidence::covering_not_crossing`].
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum SsrLocusState {
    /// The region was never walked — a bug, not data.
    NeverWalked,
    /// Walked, and no read crossed the whole tract, so none reported a length. Either no read
    /// reached the locus, or reads reached it and stopped inside it.
    NoLengthReported,
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

/// The seven selection values as a census file carries them: **one digest a field, with the
/// field's name beside it**.
///
/// **A table and not one digest over the whole**, because a mismatch has to be *named*: every
/// one of these fails the same way, silently, and "the terms differ" is not something anyone can
/// act on. **And digests rather than the values**, because [`SelectionTerms`] holds
/// `StrRepeatCriteria` and `ScanParams` — two other modules' configuration, twenty-odd scalars
/// between them — and their only use here is equality. A codec over the values would have to
/// track every field either type ever grows, and a field it failed to track would drop out of
/// the comparison silently, which is the one failure this whole check exists to prevent
/// (`arch/parameter_prepass_joint_records.md` §2.2).
///
/// **What it gives up, deliberately: a file can say *whether* it matches and not *what it was
/// built at*.** The seven values are compared and never read, so nothing is lost to a program;
/// a person wanting to know a run's seed reads it from the run.
///
/// **How a field is digested, and the direction it fails in.** The two configuration values are
/// hashed through their `Debug` rendering, so a field added anywhere inside them changes the
/// digest without anyone remembering to come back here. The cost is that a changed `Debug` impl
/// reads as changed settings — which refuses a pooling rather than allowing a wrong one, and is
/// the direction to fail in.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SelectionTermsDigest {
    /// In the order [`SelectionTerms::first_disagreement`] checks them, and named as it names
    /// them, so the sentence a mismatch produces is the same one it always was.
    fields: Vec<(&'static str, [u8; 16])>,
}

/// The seven, in the order [`SelectionTerms::first_disagreement`] checks them and under the
/// names it gives them. **One list, so the writer, the reader and the refusal cannot drift.**
pub const SELECTION_FIELDS: [&str; 7] = [
    "selection seed",
    "reference digest",
    "analysed region set",
    "repeat catalog build settings",
    "STR routing criteria",
    "generic target position count",
    "STR per-stratum cap",
];

impl SelectionTermsDigest {
    pub fn of(terms: &SelectionTerms) -> Self {
        // **Destructured without `..` on purpose.** A value added to `SelectionTerms` stops this
        // compiling rather than quietly going undigested, and a value that goes undigested lets
        // two samples asked for different loci be pooled without a word.
        let SelectionTerms {
            seed,
            reference,
            analysed_regions,
            catalog_built_under,
            ssr_criteria,
            generic_target,
            ssr_cap,
        } = terms;
        let digests = [
            digest_of(&seed.to_le_bytes()),
            digest_of(&reference.0),
            digest_of(&analysed_regions.0),
            digest_of(format!("{catalog_built_under:?}").as_bytes()),
            digest_of(format!("{ssr_criteria:?}").as_bytes()),
            digest_of(&generic_target.to_le_bytes()),
            digest_of(&(*ssr_cap as u64).to_le_bytes()),
        ];
        Self {
            fields: SELECTION_FIELDS.into_iter().zip(digests).collect(),
        }
    }

    /// Which selection value two samples first disagree on — `None` when they may be pooled.
    pub fn first_disagreement(&self, other: &Self) -> Option<&'static str> {
        self.fields
            .iter()
            .zip(&other.fields)
            .find(|(mine, theirs)| mine != theirs)
            .map(|(mine, _)| mine.0)
    }

    /// The table itself, for the codec that writes it and the one that reads it back.
    pub fn fields(&self) -> &[(&'static str, [u8; 16])] {
        &self.fields
    }

    /// A table read back from a census file.
    ///
    /// # Panics
    ///
    /// When a name is not one this build knows, or when the names do not arrive in this
    /// build's own order — either means the file was written by a build whose selection terms
    /// are not these, and comparing the two tables field by field would compare different
    /// fields.
    pub fn from_fields(fields: Vec<(&'static str, [u8; 16])>) -> Self {
        assert!(
            fields.len() == SELECTION_FIELDS.len()
                && fields
                    .iter()
                    .zip(SELECTION_FIELDS)
                    .all(|((read, _), mine)| *read == mine),
            "the selection values in this file are not the ones this build compares"
        );
        Self { fields }
    }
}

fn digest_of(bytes: &[u8]) -> [u8; 16] {
    let mut hasher = Md5::new();
    hasher.update(bytes);
    hasher.finalize().into()
}

/// The twelve values the fit refuses to pool across.
///
/// Seven say which loci were **asked for** ([`SelectionTerms`]), one says which came
/// **back** ([`CensusLociDigest`]), and five say **in what units** the evidence was written
/// down. An earlier version of this check had none of the last five, which left it able to
/// pass while two samples' rows meant different things.
#[derive(Debug, Clone, PartialEq)]
pub struct RecordingTerms {
    /// **Digested rather than held — 2026-08-14.** Their only use is a named equality, and a
    /// census file has to carry them; see [`SelectionTermsDigest`] for what that buys and what
    /// it gives up.
    pub selection: SelectionTermsDigest,
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
// What a piece of the census is, and where it sits
// ---------------------------------------------------------------------

/// Which piece of a sample's evidence — **the smallest piece anything ever asks for**.
///
/// There are exactly two kinds: one read group's ordinary-position evidence, and one read
/// group's tracts for one stratum. **That division is the estimator's own shape rather than an
/// encoding detail**: the fit finishes the ordinary positions before it reads a tract, and then
/// reads one band of strata at a time, so a run never has to hold a whole sample's evidence to
/// use any of it.
///
/// **One name for three things that must not drift**: the entries of a census file's directory,
/// the unit a scoped call lends, and the unit a counting reader measures bytes against.
///
/// **The ordering is the enumeration order**, not something a caller sorts: ordinary positions
/// before tracts, read groups by id, strata by their own key. A fit that iterates sections is
/// deterministic because this `Ord` is.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum SectionKey {
    /// One read group's ordinary-position evidence.
    Generic(ReadGroupId),
    /// One read group's tracts for one stratum.
    Ssr(ReadGroupId, Stratum),
}

impl SectionKey {
    /// The read group this section belongs to — the one thing both kinds share.
    pub fn read_group(self) -> ReadGroupId {
        match self {
            Self::Generic(group) | Self::Ssr(group, _) => group,
        }
    }

    /// The stratum, where there is one. Ordinary positions are not in a stratum.
    pub fn stratum(self) -> Option<Stratum> {
        match self {
            Self::Generic(_) => None,
            Self::Ssr(_, stratum) => Some(stratum),
        }
    }
}

/// The decoded contents of one section.
///
/// **The two halves hold different things and are read through different calls**, so this is a
/// sum rather than a shared shape: five per-base buckets cannot express a tract length and a
/// length distribution cannot express which base was substituted.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Section {
    Generic(GenericEvidence),
    Ssr(SsrEvidence),
}

impl Section {
    /// Whether this section is what the key asks for — the two must not drift apart.
    pub fn answers(&self, key: SectionKey) -> bool {
        matches!(
            (self, key),
            (Self::Generic(_), SectionKey::Generic(_)) | (Self::Ssr(_), SectionKey::Ssr(_, _))
        )
    }
}

/// One sample's ordinary-position sections, as a scoped call lends them: each read group the
/// sample holds among those asked for, with that group's section.
///
/// **The read group travels with the section because a cohort's samples need not hold the same
/// ones.** A sample sequenced in one library holds one; the fit asks for the union across the
/// cohort and each sample answers with its own, which is what it did when the maps were public.
///
/// **A borrow with no owner of its own.** The whole point of the scoped calls is that this
/// cannot outlive the call that produced it, so it is a `Vec` of borrows and never a `Vec` of
/// sections.
pub type SampleGenericSections<'a> = Vec<(ReadGroupId, &'a GenericEvidence)>;

/// Where a section sits in a file: an offset and a length.
///
/// **A pair and not a range**, because it is written into a directory and read back as two
/// numbers, and because a length says what to allocate where an end does not.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ByteExtent {
    offset: u64,
    len: u64,
}

impl ByteExtent {
    pub fn new(offset: u64, len: u64) -> Self {
        Self { offset, len }
    }

    pub fn offset(self) -> u64 {
        self.offset
    }

    pub fn len(self) -> u64 {
        self.len
    }

    pub fn is_empty(self) -> bool {
        self.len == 0
    }

    /// Whether two sections overlap — **which they never may**, since a call reading one seeks
    /// to its own offset and takes its own length.
    pub fn overlaps(self, other: Self) -> bool {
        self.offset < other.offset + other.len && other.offset < self.offset + self.len
    }
}

/// Where a sample's sections are.
///
/// **Two states and no more**: a genome walk produces the first, opening a census file produces
/// the second, and a caller cannot tell which it has — that is what lets the fit be one code
/// path over a cohort held in memory and a cohort read from disk.
///
/// **Only the resident state exists today.** The file, its directory and the reader that seeks
/// into it are the next unit of work; nothing here changes when they arrive, which is why this
/// type is written now.
///
/// **This type's visibility is not what keeps a section from being retained.** What does that
/// is the scoped access on the value that owns one: a section is lent for the length of a call
/// and there is no field a caller could keep it in.
#[derive(Debug, Clone, PartialEq)]
pub enum Sections {
    /// Built by a walk: every section decoded and held. This is the run that never writes a
    /// file, at fifty samples or at one.
    Resident(BTreeMap<SectionKey, Section>),
    /// Opened and checked, the directory read, **nothing else decoded**. A call fills the
    /// sections it was asked for, hands them over, and drops them when it returns — there is no
    /// field here one could be kept in.
    ///
    /// **A path and not an open reader — 2026-08-14, and it is a deviation worth its reason.**
    /// The architecture (§1.1a) holds a `Box<dyn ReadSeek>`, which would make a cohort hold one
    /// open file a sample: this caller commits to cohorts of **several thousand**
    /// (`design_principles.md` §0), and a thousand open descriptors sits at or above the default
    /// soft limit on both platforms this runs on — before the pileups are counted. A path opens
    /// the file once a *call* instead, seeks to each section inside it, and closes it again; the
    /// contract §1.1a states — one read a section, decoded from a slice, nothing retained
    /// between calls — is unchanged, and a value that names a file rather than holding one open
    /// can also be cloned and compared, which the fit's own tests and examples need.
    Backed {
        path: PathBuf,
        directory: BTreeMap<SectionKey, ByteExtent>,
    },
}

impl Sections {
    /// The sections a walk built, checked against their own keys.
    ///
    /// # Panics
    ///
    /// When a key and its section disagree about which half of the census they are — an
    /// ordinary-position key holding tracts indexes one thing by another's positions, and
    /// nothing downstream would notice.
    pub fn resident(sections: BTreeMap<SectionKey, Section>) -> Self {
        assert!(
            sections.iter().all(|(key, section)| section.answers(*key)),
            "a section is filed under a key of the other kind"
        );
        Self::Resident(sections)
    }

    /// A census file, opened and its directory read — **and nothing else decoded**.
    pub fn backed(path: PathBuf, directory: BTreeMap<SectionKey, ByteExtent>) -> Self {
        Self::Backed { path, directory }
    }

    /// Every section this sample holds, in enumeration order.
    ///
    /// **A file-backed value answers this from its directory**, which is why the directory is at
    /// the head of the file: knowing what is there costs a few hundred bytes rather than the
    /// whole of it.
    pub fn keys(&self) -> Box<dyn Iterator<Item = SectionKey> + '_> {
        match self {
            Self::Resident(sections) => Box::new(sections.keys().copied()),
            Self::Backed { directory, .. } => Box::new(directory.keys().copied()),
        }
    }

    /// Whether this sample holds the section named — answered without decoding one.
    pub fn holds(&self, key: SectionKey) -> bool {
        match self {
            Self::Resident(sections) => sections.contains_key(&key),
            Self::Backed { directory, .. } => directory.contains_key(&key),
        }
    }

    /// Fill `into` with the sections named by `keys`, **in that order**, for a caller about to
    /// borrow them.
    ///
    /// **Resident leaves it empty**: what it holds is already decoded, and copying it would be
    /// the memory this design exists to avoid. **Backed opens the file once, seeks to each
    /// section, reads exactly its own bytes and decodes from the slice** — which is the
    /// contract, not a preference: reading fields through the file handle would pay an
    /// indirection per byte and, unbuffered, a syscall per field.
    ///
    /// # Errors
    ///
    /// [`CensusError::Io`] when the file will not open or read, [`CensusError::Malformed`] when
    /// a section's bytes are not a section or a key is not in the directory.
    pub(super) fn fill(
        &self,
        keys: &[SectionKey],
        into: &mut Vec<Section>,
    ) -> Result<(), CensusError> {
        into.clear();
        let Self::Backed { path, directory } = self else {
            return Ok(());
        };
        let mut file = std::fs::File::open(path)?;
        let mut buffer = Vec::new();
        for key in keys {
            let extent = directory.get(key).ok_or(CensusError::Malformed)?;
            file.seek(std::io::SeekFrom::Start(extent.offset()))?;
            let len = usize::try_from(extent.len()).map_err(|_| CensusError::Malformed)?;
            buffer.resize(len, 0);
            file.read_exact(&mut buffer)?;
            // Counted where the bytes actually leave the file, so what spec §7.15 asserts is the
            // read itself and not a promise about it.
            super::census_file::count_bytes_read(extent.len());
            into.push(super::census_file::decode_section(*key, &buffer)?);
        }
        Ok(())
    }

    /// The sections named by `keys`, in that order — from what this value holds, or from what
    /// [`fill`](Self::fill) just decoded.
    ///
    /// **Neither is a section this call could hand back**: the borrows live as long as `filled`
    /// does, and `filled` is a local of the scoped call above.
    ///
    /// # Panics
    ///
    /// When a key is one this sample does not hold, which [`fill`](Self::fill) has already
    /// refused for a file — so reaching it means a resident caller asked for a section that was
    /// never recorded, and that is a selection disagreement rather than a file.
    pub(super) fn borrow<'a>(
        &'a self,
        keys: &[SectionKey],
        filled: &'a [Section],
    ) -> Vec<&'a Section> {
        match self {
            Self::Resident(sections) => keys
                .iter()
                .map(|key| {
                    sections
                        .get(key)
                        .unwrap_or_else(|| panic!("this census holds no section for {key:?}"))
                })
                .collect(),
            Self::Backed { .. } => filled.iter().collect(),
        }
    }

    /// Every section with its key, in enumeration order.
    ///
    /// **Resident only, and inside this module only.** A file-backed value cannot answer this
    /// without decoding the whole file, which is the thing it exists not to do; what it will
    /// answer is [`keys`](Self::keys) and one section at a time.
    pub(super) fn iter(&self) -> impl Iterator<Item = (SectionKey, &Section)> + '_ {
        match self {
            Self::Resident(sections) => sections.iter().map(|(key, section)| (*key, section)),
            Self::Backed { .. } => {
                panic!("a file-backed census is read one section at a time, never all at once")
            }
        }
    }

    pub fn len(&self) -> usize {
        match self {
            Self::Resident(sections) => sections.len(),
            Self::Backed { directory, .. } => directory.len(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

// ---------------------------------------------------------------------
// The whole input to the fit
// ---------------------------------------------------------------------

/// One sample's evidence at the kept loci, plus the values the fit checks before pooling.
///
/// **The sections are not public, and that is the whole point of the type — 2026-08-14.** A
/// caller that could reach the map could keep every section it ever looked at, at which point a
/// run reading a census file has quietly reassembled the file in memory — the outcome the
/// by-section design exists to prevent, arrived at without anyone deciding to
/// (`parameter_prepass_joint_records.md` §6.2). So a section is **lent for the length of a
/// call** by [`with_generic`](Self::with_generic) and [`with_strata`](Self::with_strata), and
/// there is no field a caller could put one in.
///
/// **One type for both kinds of run.** A walk fills every section and holds it; a file-backed
/// value will fill one when it is asked for and drop it when the call returns. That difference
/// lives in [`Sections`], so the parameters fit is written once.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleCensusEvidence {
    pub sample: String,
    pub terms: RecordingTerms,
    sections: Sections,
}

impl SampleCensusEvidence {
    /// A sample whose sections a walk built and which are all in memory.
    pub fn resident(
        sample: String,
        terms: RecordingTerms,
        sections: BTreeMap<SectionKey, Section>,
    ) -> Self {
        Self {
            sample,
            terms,
            sections: Sections::resident(sections),
        }
    }

    /// A sample whose sections are a file, opened and its directory read — **and nothing else
    /// decoded**. What a scoped call asks for is filled then, and dropped when it returns.
    pub fn backed(
        sample: String,
        terms: RecordingTerms,
        path: PathBuf,
        directory: BTreeMap<SectionKey, ByteExtent>,
    ) -> Self {
        Self {
            sample,
            terms,
            sections: Sections::backed(path, directory),
        }
    }

    /// Every read group this sample recorded, in the order the fit visits them.
    pub fn read_groups(&self) -> Vec<ReadGroupId> {
        self.sections
            .keys()
            .map(|key| key.read_group())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    /// Every stratum this sample holds tracts for, in the stratum's own order.
    pub fn strata(&self) -> Vec<Stratum> {
        self.sections
            .keys()
            .filter_map(|key| key.stratum())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect()
    }

    /// This sample's ordinary-position evidence for a band of read groups, **for the length of
    /// the call**.
    ///
    /// **A band and not one read group**, for the reason [`with_strata`](Self::with_strata)
    /// takes a band of strata: a position's depth is the sum of the sample's read groups'
    /// depths, so a call lending one group at a time would leave the fit holding a partial sum
    /// for every kept position while it asked for the next.
    ///
    /// # Errors
    ///
    /// [`CensusError`] where this sample is a file and the file will not answer.
    ///
    /// # Panics
    ///
    /// When a *resident* sample holds no section for one of the groups — a stratum or a group
    /// with no reads is still a section holding zero, so an absent one means the caller and the
    /// census were built from different selections.
    pub fn with_generic<R>(
        &mut self,
        groups: &[ReadGroupId],
        f: impl FnOnce(&[&GenericEvidence]) -> R,
    ) -> Result<R, CensusError> {
        let keys: Vec<SectionKey> = groups.iter().map(|g| SectionKey::Generic(*g)).collect();
        let mut filled = Vec::new();
        self.sections.fill(&keys, &mut filled)?;
        let sections: Vec<&GenericEvidence> = self
            .sections
            .borrow(&keys, &filled)
            .into_iter()
            .map(|section| match section {
                Section::Generic(records) => records,
                Section::Ssr(_) => {
                    panic!("a section filed under an ordinary-position key holds tracts")
                }
            })
            .collect();
        Ok(f(&sections))
    }

    /// This sample's tracts for one read group over a band of strata, **for the length of the
    /// call**.
    ///
    /// # Errors
    ///
    /// As [`with_generic`](Self::with_generic).
    ///
    /// # Panics
    ///
    /// When a resident sample holds no section for one of the strata — see
    /// [`with_generic`](Self::with_generic).
    pub fn with_strata<R>(
        &mut self,
        group: ReadGroupId,
        strata: &[Stratum],
        f: impl FnOnce(&[&SsrEvidence]) -> R,
    ) -> Result<R, CensusError> {
        let keys: Vec<SectionKey> = strata
            .iter()
            .map(|stratum| SectionKey::Ssr(group, *stratum))
            .collect();
        let mut filled = Vec::new();
        self.sections.fill(&keys, &mut filled)?;
        let sections: Vec<&SsrEvidence> = self
            .sections
            .borrow(&keys, &filled)
            .into_iter()
            .map(|section| match section {
                Section::Ssr(records) => records,
                Section::Generic(_) => {
                    panic!("a section filed under a tract key holds ordinary positions")
                }
            })
            .collect();
        Ok(f(&sections))
    }

    /// This sample's depth at the `index`-th kept generic position, summed over its read
    /// groups.
    ///
    /// **Summing is where the two grains meet.** The cross-sample work wants the *sample's*
    /// counts at a position; the error rate is chemistry and wants the *read group's*. The
    /// record is kept at the finer grain because summing is addition of raw counts at one
    /// place and the reverse is not recoverable.
    pub fn allele_counts(&self, index: usize) -> [u32; 5] {
        let mut counts = [0_u32; 5];
        for (_, records) in self.generic_sections() {
            for observation in records.at(index).1 {
                counts[observation.allele.code() as usize] += u32::from(observation.reads);
            }
        }
        counts
    }

    /// Every ordinary-position section, with the read group it belongs to.
    ///
    /// **Inside this module only, and it is not the door the fit comes in through.** What keeps
    /// a section from being retained is the scoped calls above; this is the primitive they and
    /// [`CohortCensusEvidence`] are built on, because a scoped call cannot be composed over a
    /// list of samples whose length is not known while the code is being written.
    pub(super) fn generic_sections(&self) -> Vec<(ReadGroupId, &GenericEvidence)> {
        self.sections
            .iter()
            .filter_map(|(key, section)| match (key, section) {
                (SectionKey::Generic(group), Section::Generic(records)) => Some((group, records)),
                _ => None,
            })
            .collect()
    }

    /// The keys this sample holds among the ordinary-position sections asked for.
    pub(super) fn generic_keys(&self, groups: &[ReadGroupId]) -> Vec<SectionKey> {
        groups
            .iter()
            .map(|group| SectionKey::Generic(*group))
            .filter(|key| self.sections.holds(*key))
            .collect()
    }

    /// The keys this sample holds among the tract sections in a band of strata, over every read
    /// group it recorded.
    pub(super) fn ssr_keys(&self, strata: &[Stratum]) -> Vec<SectionKey> {
        self.sections
            .keys()
            .filter(|key| key.stratum().is_some_and(|s| strata.contains(&s)))
            .collect()
    }

    /// Fill `into` with the sections named, for a caller about to borrow them — see
    /// [`Sections::fill`].
    pub(super) fn fill(
        &self,
        keys: &[SectionKey],
        into: &mut Vec<Section>,
    ) -> Result<(), CensusError> {
        self.sections.fill(keys, into)
    }

    /// This sample's ordinary-position sections for `keys`, borrowed from what it holds or from
    /// what [`fill`](Self::fill) decoded.
    pub(super) fn borrow_generic<'a>(
        &'a self,
        keys: &[SectionKey],
        filled: &'a [Section],
    ) -> SampleGenericSections<'a> {
        keys.iter()
            .zip(self.sections.borrow(keys, filled))
            .filter_map(|(key, section)| match (key, section) {
                (SectionKey::Generic(group), Section::Generic(records)) => Some((*group, records)),
                _ => None,
            })
            .collect()
    }

    /// The same on the tract half.
    pub(super) fn borrow_ssr<'a>(
        &'a self,
        keys: &[SectionKey],
        filled: &'a [Section],
    ) -> SampleTractSections<'a> {
        keys.iter()
            .zip(self.sections.borrow(keys, filled))
            .filter_map(|(key, section)| match (key, section) {
                (SectionKey::Ssr(group, stratum), Section::Ssr(records)) => {
                    Some((*group, *stratum, records))
                }
                _ => None,
            })
            .collect()
    }

    /// Every repeat-tract section, with the read group and stratum it belongs to. See
    /// [`generic_sections`](Self::generic_sections) for why this exists.
    pub(super) fn ssr_sections(&self) -> Vec<(ReadGroupId, Stratum, &SsrEvidence)> {
        self.sections
            .iter()
            .filter_map(|(key, section)| match (key, section) {
                (SectionKey::Ssr(group, stratum), Section::Ssr(records)) => {
                    Some((group, stratum, records))
                }
                _ => None,
            })
            .collect()
    }
}

/// One sample's repeat-tract sections, as a scoped call lends them: each (read group, stratum)
/// the sample holds among the strata asked for.
///
/// **The stratum travels with the section** because a call lends a *band* of them — the fit
/// borrows a thin stratum together with its neighbours (spec §6.2) — so the closure has to be
/// able to tell which is which.
pub type SampleTractSections<'a> = Vec<(ReadGroupId, Stratum, &'a SsrEvidence)>;

/// Two samples that did not record the same thing, and the value they first differ on.
///
/// **Naming the value is the whole point** (spec §5): every value in [`RecordingTerms`] fails
/// the same way, silently, and only its name says what to fix.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TermsDisagreement {
    pub first: String,
    pub second: String,
    /// What the two disagree on, in the words [`RecordingTerms::first_disagreement`] uses.
    pub field: &'static str,
}

impl std::fmt::Display for TermsDisagreement {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "samples {} and {} disagree on {}; they did not record the same thing",
            self.first, self.second, self.field
        )
    }
}

/// Every sample's census, however each one is held.
///
/// **This is what the parameters fit reads, and it reads it a band at a time.** A tract's length
/// frequencies are fitted from every sample with reads there, so the unit lent is one band of
/// strata across the *whole cohort* rather than one sample's stratum (spec §6.2) — which is why
/// the scoped calls live here and not only on the sample.
///
/// **The twelve recording terms are checked before a single section is read.** They say which
/// loci were asked for, which came back, and in what units the evidence was written down; two
/// samples that disagree on any of them hold rows that mean different things, and every one of
/// them fails silently. So [`new`](Self::new) refuses the cohort at the door, where the refusal
/// costs nothing, rather than after a pass over two million positions.
#[derive(Debug, Clone, PartialEq)]
pub struct CohortCensusEvidence {
    samples: Vec<SampleCensusEvidence>,
    /// Every read group any sample recorded, in the order the fit visits them. **The union and
    /// not one sample's**: a cohort's samples need not have been sequenced the same way.
    read_groups: Vec<ReadGroupId>,
    /// Every stratum any sample holds tracts for. These do agree across samples — they come
    /// from the selection, which the terms above have already been checked on.
    strata: Vec<Stratum>,
}

impl CohortCensusEvidence {
    /// Adopt each sample's census, refusing the cohort if two of them did not record the same
    /// thing.
    ///
    /// # Errors
    ///
    /// [`TermsDisagreement`], naming the first sample that disagrees with the first and the
    /// value they differ on. **Checked before any section is decoded**, which is what makes the
    /// refusal free.
    pub fn new(samples: Vec<SampleCensusEvidence>) -> Result<Self, TermsDisagreement> {
        if let Some(first) = samples.first() {
            for other in &samples[1..] {
                if let Some(field) = first.terms.first_disagreement(&other.terms) {
                    return Err(TermsDisagreement {
                        first: first.sample.clone(),
                        second: other.sample.clone(),
                        field,
                    });
                }
            }
        }
        let read_groups = samples
            .iter()
            .flat_map(|sample| sample.read_groups())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        let strata = samples
            .iter()
            .flat_map(|sample| sample.strata())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        Ok(Self {
            samples,
            read_groups,
            strata,
        })
    }

    pub fn len(&self) -> usize {
        self.samples.len()
    }

    pub fn is_empty(&self) -> bool {
        self.samples.is_empty()
    }

    /// Each sample's name, in the order every scoped call indexes by.
    pub fn sample_names(&self) -> impl Iterator<Item = &str> + '_ {
        self.samples.iter().map(|sample| sample.sample.as_str())
    }

    /// The terms every sample in this cohort recorded under — one set, since
    /// [`new`](Self::new) has already refused any cohort where they differ.
    pub fn terms(&self) -> Option<&RecordingTerms> {
        self.samples.first().map(|sample| &sample.terms)
    }

    pub fn read_groups(&self) -> &[ReadGroupId] {
        &self.read_groups
    }

    pub fn strata(&self) -> &[Stratum] {
        &self.strata
    }

    /// Every sample's ordinary-position evidence for a band of read groups, **for the length of
    /// the call**. Entry `s` is sample `s`'s own sections among those asked for — a sample
    /// sequenced in one library answers with one where another answers with three.
    pub fn with_generic<R>(
        &mut self,
        groups: &[ReadGroupId],
        f: impl FnOnce(&[SampleGenericSections<'_>]) -> R,
    ) -> Result<R, CensusError> {
        let keys: Vec<Vec<SectionKey>> = self
            .samples
            .iter()
            .map(|sample| sample.generic_keys(groups))
            .collect();
        // **Decoded first, borrowed second**, and the buffers are locals of this call: when it
        // returns they go, and a file-backed sample has nothing left decoded.
        let mut filled: Vec<Vec<Section>> = (0..self.samples.len()).map(|_| Vec::new()).collect();
        for ((sample, keys), into) in self.samples.iter().zip(&keys).zip(&mut filled) {
            sample.fill(keys, into)?;
        }
        let lent: Vec<SampleGenericSections<'_>> = self
            .samples
            .iter()
            .zip(&keys)
            .zip(&filled)
            .map(|((sample, keys), filled)| sample.borrow_generic(keys, filled))
            .collect();
        Ok(f(&lent))
    }

    /// Every sample's tracts for a band of strata, **for the length of the call**.
    ///
    /// **A band and not one stratum**: 68 of tomato's 141 strata hold fewer than a hundred
    /// tracts and are fitted by borrowing from their neighbouring repeat counts
    /// (`parameter_prepass_joint_loci.md` §3.6), so the fit sometimes wants a stratum and its
    /// neighbours together.
    pub fn with_strata<R>(
        &mut self,
        strata: &[Stratum],
        f: impl FnOnce(&[SampleTractSections<'_>]) -> R,
    ) -> Result<R, CensusError> {
        let keys: Vec<Vec<SectionKey>> = self
            .samples
            .iter()
            .map(|sample| sample.ssr_keys(strata))
            .collect();
        let mut filled: Vec<Vec<Section>> = (0..self.samples.len()).map(|_| Vec::new()).collect();
        for ((sample, keys), into) in self.samples.iter().zip(&keys).zip(&mut filled) {
            sample.fill(keys, into)?;
        }
        let lent: Vec<SampleTractSections<'_>> = self
            .samples
            .iter()
            .zip(&keys)
            .zip(&filled)
            .map(|((sample, keys), filled)| sample.borrow_ssr(keys, filled))
            .collect();
        Ok(f(&lent))
    }

    /// The samples themselves, for a caller that has to write them down.
    ///
    /// **This is not a way to reach a section**: a `SampleCensusEvidence` keeps its own private,
    /// and the scoped calls are still the only door.
    pub fn samples(&self) -> &[SampleCensusEvidence] {
        &self.samples
    }

    /// The samples back, for a caller that has finished with the cohort.
    pub fn into_samples(self) -> Vec<SampleCensusEvidence> {
        self.samples
    }
}

/// Why a census could not be written or read.
///
/// **The guard threshold is not in here.** A locus over one non-whole-repeat read in ten is
/// well-formed data the parameters fit should decline to fit (spec §3.3); encoding it as a write
/// failure would make a property of the sample look like a property of the file.
#[non_exhaustive]
#[derive(Debug, thiserror::Error)]
pub enum CensusError {
    /// The stream is not a census, or is a version this build does not read.
    #[error("the stream is not a census, or is of a version this build does not read")]
    Malformed,
    #[error("i/o while reading or writing a census")]
    Io(#[from] std::io::Error),
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
    /// Which stratum each of those loci is in, in the same order. **The writer carries it
    /// because the tract half is stored one stratum at a time**, and the locus stream says which
    /// tract a read reached, never which stratum it belongs to.
    ssr_stratum_of: Vec<Stratum>,
    /// Where each of those loci sits **within its own stratum**, in the same order — the index
    /// every entry of an [`SsrEvidence`] carries, since a section holds one stratum's tracts and
    /// a genome-wide index would name a tract it does not have.
    ssr_index_in_stratum: Vec<u32>,
    /// How many tracts each stratum holds, which is the length a section of it is built at
    /// whether or not a read ever reaches one of them.
    ssr_stratum_size: BTreeMap<Stratum, usize>,
    edges: DepthBinEdges,
    generic: BTreeMap<ReadGroupId, GenericEvidence>,
    /// Per read group, then per stratum. **Two levels rather than one key**, because the
    /// walk asks twice for "every read group seen so far" and the outer map answers it.
    ssr: BTreeMap<ReadGroupId, BTreeMap<Stratum, SsrEvidence>>,
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
        // **The stratum travels beside the locus from here on**, because the locus stream the
        // walk hands over names a tract by where it sits and never by which stratum it is in,
        // and two of a tract's counts are kept per stratum.
        let mut ssr_loci: Vec<(GenomeRegion, Stratum)> = loci
            .ssr()
            .iter_sorted()
            .into_iter()
            .flat_map(|((period, reference_repeats), segments)| {
                segments.iter().map(move |segment| {
                    (
                        segment,
                        Stratum {
                            period,
                            reference_repeats,
                        },
                    )
                })
            })
            .filter_map(|(segment, stratum)| {
                contig_of(segment.chrom()).map(|contig| {
                    (
                        GenomeRegion {
                            contig,
                            start: Position(segment.start()),
                            end: Position(segment.end()),
                        },
                        stratum,
                    )
                })
            })
            .collect();
        ssr_loci.sort_unstable_by_key(|(r, _)| (r.contig.get(), r.start.get(), r.end.get()));
        let ssr_stratum_of: Vec<Stratum> = ssr_loci.iter().map(|(_, stratum)| *stratum).collect();
        // **Numbered within its stratum in genome order**, so that a section's tracts run in the
        // same order the whole selection does. The fit sums over a stratum's tracts, and a sum of
        // floats is not associative — an order that depended on how the strata were interleaved
        // would move the answer's last digits for no reason anybody could name.
        let mut ssr_stratum_size: BTreeMap<Stratum, usize> = BTreeMap::new();
        let mut ssr_index_in_stratum: Vec<u32> = Vec::with_capacity(ssr_stratum_of.len());
        for stratum in &ssr_stratum_of {
            let so_far = ssr_stratum_size.entry(*stratum).or_insert(0);
            ssr_index_in_stratum.push(*so_far as u32);
            *so_far += 1;
        }
        let ssr_loci = ssr_loci.into_iter().map(|(region, _)| region).collect();
        Self {
            generic_loci: loci.generic().to_vec(),
            ssr: BTreeMap::new(),
            ssr_loci,
            ssr_stratum_of,
            ssr_index_in_stratum,
            ssr_stratum_size,
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
        let stratum = self.ssr_stratum_of[index];
        // The index every entry below carries: this tract's place **in its own stratum**, since
        // that is the whole of what a section holds.
        let local = self.ssr_index_in_stratum[index] as usize;
        let tracts_here = self.ssr_stratum_size[&stratum];
        // **How many of this read group's reads at this locus have been numbered already.**
        // The read a difference sits on is numbered within the *locus*, not within the
        // observation it arrived in: the walk folds reads carrying identical bases into one
        // observation, so a counter that restarted at each observation gave two reads that
        // differ from the reference in different places the same number, and a consumer
        // grouping by it saw one bad read where the truth is an allele on two reads.
        //
        // It advances over the reads that were **compared** — the unslipped ones that reach
        // the difference loop below — because those are the only reads that can carry a
        // difference, and every number spent elsewhere is headroom lost from a byte.
        let mut numbered: BTreeMap<ReadGroupId, u32> = BTreeMap::new();

        for observation in &locus.observations {
            let records = self
                .ssr
                .entry(observation.read_group)
                .or_default()
                .entry(stratum)
                .or_insert_with(|| SsrEvidence::never_walked(tracts_here));
            records.walked.set(local);
            let reads = u16::try_from(observation.num_obs).unwrap_or(u16::MAX);
            if observation.read_witness != crate::ng::locus_generation::ReadWitness::Complete {
                // A read that covered the tract without crossing it reports no length; it is a
                // lower bound. **Counted into this stratum's total rather than at the locus**
                // (spec §3): the fit summed the per-locus version the moment it read it, and
                // the loss runs along repeat count, which is what a stratum is.
                records.covering_not_crossing += u64::from(reads);
                continue;
            }
            let difference = observation.bases.len() as i64 - reference_length;
            if difference % period == 0 {
                records.offsets[local].add((difference / period) as i32, reads);
            } else {
                records.guard.push(GuardObservation {
                    locus: local as u32,
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
            // **Per stratum, for the same reason**: it is the denominator of a rate fitted per
            // stratum, and per locus it was 1.85 MB a read group on tomato that nothing read.
            records.bases_compared += u64::from(reads) * observation.bases.len() as u64;
            // **One entry per read, not one per distinct sequence.** The walk has already
            // folded reads carrying the same bases into a single observation with a count, and
            // collapsing them here too would lose the fact this channel exists for: that two
            // reads carried the *same* interruption, which is what makes an interruption an
            // allele rather than an error.
            let numbered_here = numbered.entry(observation.read_group).or_insert(0);
            let first_read = *numbered_here;
            *numbered_here += u32::from(reads);
            for (offset, (read_base, reference_base)) in observation
                .bases
                .iter()
                .zip(locus.reference_bases.iter())
                .enumerate()
            {
                if read_base.eq_ignore_ascii_case(reference_base) {
                    continue;
                }
                for copy in 0..u32::from(reads) {
                    records.differences.push(TractDifference {
                        locus: local as u32,
                        read: u16::try_from(first_read + copy).unwrap_or(u16::MAX),
                        offset: i16::try_from(offset).unwrap_or(i16::MAX),
                        base: ObservedAllele::of_base(*read_base),
                    });
                }
            }
        }
        // Reads that reached the locus and produced no observation at all. **Every read group
        // seen so far learns the walk reached this tract**, and the count itself is charged
        // once — to the lowest-numbered group, since the reads cannot be attributed to one and
        // the fit adds every group's total into the stratum's anyway.
        if locus.reads_without_observation > 0 {
            for strata in self.ssr.values_mut() {
                strata
                    .entry(stratum)
                    .or_insert_with(|| SsrEvidence::never_walked(tracts_here))
                    .walked
                    .set(local);
            }
            if let Some(strata) = self.ssr.values_mut().next() {
                strata
                    .entry(stratum)
                    .or_insert_with(|| SsrEvidence::never_walked(tracts_here))
                    .covering_not_crossing += u64::from(locus.reads_without_observation);
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

        // **The same sections in every sample, and that is what a census is.** Every declared
        // read group gets an ordinary-position section and one section per stratum the selection
        // holds, whether or not a read ever reached one of them — a stratum with no tracts read
        // is a section holding zero, and an absent section would say instead that this sample was
        // never asked. The fit enumerates sections; two samples enumerating different ones would
        // be two samples answering different questions.
        let mut sections: BTreeMap<SectionKey, Section> = BTreeMap::new();
        for (group, records) in std::mem::take(&mut self.generic) {
            sections.insert(SectionKey::Generic(group), Section::Generic(records));
        }
        let mut tracts = std::mem::take(&mut self.ssr);
        for group in &self.read_groups {
            let mut theirs = tracts.remove(group).unwrap_or_default();
            for (stratum, size) in &self.ssr_stratum_size {
                let records = theirs
                    .remove(stratum)
                    .unwrap_or_else(|| SsrEvidence::never_walked(*size));
                sections.insert(SectionKey::Ssr(*group, *stratum), Section::Ssr(records));
            }
        }
        // **The declared groups are all the groups there are.** The generic half already drops
        // an observation whose read group the sample did not declare, so tract evidence from one
        // would be half a read group — recorded at the tracts and absent at the positions its
        // own error rate is fitted from.
        assert!(
            tracts.is_empty(),
            "read group {} put reads on a repeat tract without being declared by the sample",
            tracts.keys().next().expect("a group").get()
        );

        SampleCensusEvidence::resident(
            self.sample,
            RecordingTerms {
                selection: SelectionTermsDigest::of(&self.terms),
                kept_loci: self.digester.finish(),
                ssr_stratum_counts: self.stratum_counts,
                read_cap: self.read_cap,
                depth_ladder: DepthLadderDigest::of(&self.edges),
                depth_cap: self.depth_cap,
            },
            sections,
        )
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
        let records = SampleCensusEvidence::resident(
            "s".to_string(),
            terms(),
            BTreeMap::from([
                (SectionKey::Generic(ReadGroupId(0)), Section::Generic(one)),
                (SectionKey::Generic(ReadGroupId(1)), Section::Generic(two)),
            ]),
        );
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

    /// A bit a locus, set independently at every offset within a byte — the mutation this
    /// catches is a mask that also clears its neighbours, which every other test in the module
    /// passes through because they walk one locus at a time.
    #[test]
    fn one_locus_marked_walked_leaves_every_other_locus_alone() {
        let mut walked = WalkedBits::none_of(20);
        assert!((0..20).all(|locus| !walked.get(locus)));
        assert_eq!(
            walked.as_bytes().len(),
            3,
            "twenty loci is three bytes, where a bool a locus is twenty"
        );

        for locus in (0..20).step_by(3) {
            walked.set(locus);
        }
        for locus in 0..20 {
            assert_eq!(walked.get(locus), locus % 3 == 0, "locus {locus}");
        }
    }

    /// **Spec §7.4, and the answer it left open.** The fourth state — reads reached the tract
    /// and crossed none of it — was to be planted deliberately, so that the test "either shows
    /// a lower bound surviving or shows that it does not". It does not survive **per locus**:
    /// the count moved to the stratum on 2026-08-14, so a locus whose reads all stopped inside
    /// the tract reads the same as one no read reached. Both reported no length, which is the
    /// only thing the fit does with either. §7.4's other half is the test below it, where the
    /// state does survive — at the grain it is acted on.
    #[test]
    fn the_three_states_at_an_str_locus_are_distinguishable() {
        let mut records = SsrEvidence::never_walked(3);
        assert_eq!(records.state(0), SsrLocusState::NeverWalked);

        records.walked.set(1);
        assert_eq!(records.state(1), SsrLocusState::NoLengthReported);

        records.walked.set(2);
        records.offsets[2].add(0, 5);
        assert_eq!(records.state(2), SsrLocusState::Crossed);
    }

    /// The lower bound survives **at the stratum**, which is where the censoring runs: a tract
    /// longer than the reads is never crossed in any sample at any depth, so the loss follows
    /// repeat count, and repeat count is half of what a stratum is.
    #[test]
    fn a_read_that_reached_a_tract_without_crossing_it_is_counted_against_its_stratum() {
        let mut writer = ssr_writer();
        let mut covering = observation(b"ATATAT", 0, 3);
        covering.read_witness =
            ReadWitness::from_left(6, crate::ng::locus_generation::LocusLen::from_positions(12))
                .expect("a read reaching six of the tract's twelve bases");
        writer.add_locus(&ssr_locus(vec![
            covering,
            observation(b"ATATATATATAT", 0, 2),
        ]));
        let evidence = writer.finish();
        let ssr = ssr_of(&evidence, 0, AT_SIX_REPEATS);

        assert_eq!(
            ssr.covering_not_crossing(),
            3,
            "the three reads that stopped inside the tract are this stratum's censored ones"
        );
        assert_eq!(
            ssr.state(0),
            SsrLocusState::Crossed,
            "the locus itself reports the length the two crossing reads showed"
        );
    }

    /// **The count is only worth moving if it lands in the right stratum**, which is the one
    /// thing a per-locus count could not get wrong. Two tracts, two strata, one read group:
    /// each stratum is charged its own censored reads and its own compared bases.
    #[test]
    fn two_strata_are_charged_separately_for_their_censored_reads_and_compared_bases() {
        let mut writer = two_stratum_writer();

        let mut short = observation(b"ATATAT", 0, 3);
        short.read_witness =
            ReadWitness::from_left(6, crate::ng::locus_generation::LocusLen::from_positions(12))
                .expect("a read reaching six of the tract's twelve bases");
        writer.add_locus(&ssr_locus(vec![short, observation(b"ATATATATATAT", 0, 2)]));

        let mut stopped = observation(b"ATGATG", 0, 5);
        stopped.read_witness =
            ReadWitness::from_left(6, crate::ng::locus_generation::LocusLen::from_positions(12))
                .expect("a read reaching six of the tract's twelve bases");
        writer.add_locus(&atg_locus(vec![
            stopped,
            observation(b"ATGATGATGATG", 0, 4),
        ]));

        let evidence = writer.finish();
        // **Two strata, two sections**, and each one's counts are its own — since 2026-08-14 a
        // section *is* one stratum, so what a wrong stratum would do is put the count in the
        // other section rather than under the other key of one map.
        let at = ssr_of(&evidence, 0, AT_SIX_REPEATS);
        let atg = ssr_of(&evidence, 0, ATG_FOUR_REPEATS);

        assert_eq!(at.covering_not_crossing(), 3);
        assert_eq!(atg.covering_not_crossing(), 5);
        assert_eq!(
            at.bases_compared(),
            2 * 12,
            "two crossing reads over a twelve-base tract"
        );
        assert_eq!(
            atg.bases_compared(),
            4 * 12,
            "four of them at the other tract, and neither total leaks into the other"
        );
    }

    /// **The index a tract's entries carry is its place in its own stratum, and nothing else.**
    /// Two twelve-base tracts of different motifs: the second is the second locus of the
    /// selection and the *first* of its own stratum, so a writer that kept the genome-wide
    /// index would write a guard entry and a mismatch against tract 1 of a one-tract section.
    #[test]
    fn a_tract_is_numbered_within_its_stratum_and_not_across_the_selection() {
        let mut writer = two_stratum_writer();
        // The `AT` tract, clean and at the reference length.
        writer.add_locus(&ssr_locus(vec![observation(b"ATATATATATAT", 0, 2)]));
        // The `ATG` tract: one read carrying an interruption, and one whose length differs by
        // a base — not a whole copy, so the guard catches it.
        writer.add_locus(&atg_locus(vec![
            observation(b"ATGATGCTGATG", 0, 1),
            observation(b"ATGATGATGATGA", 0, 1),
        ]));
        let evidence = writer.finish();

        let at = ssr_of(&evidence, 0, AT_SIX_REPEATS);
        let atg = ssr_of(&evidence, 0, ATG_FOUR_REPEATS);
        assert_eq!((at.len(), atg.len()), (1, 1), "one tract in each stratum");
        assert_eq!(
            atg.differences()
                .iter()
                .map(|d| d.locus)
                .collect::<Vec<_>>(),
            vec![0],
            "the second locus of the selection is the first of its own stratum"
        );
        assert_eq!(
            atg.guard().iter().map(|g| g.locus).collect::<Vec<_>>(),
            vec![0],
            "and so is the guard's"
        );
        assert!(at.differences().is_empty() && at.guard().is_empty());
        assert_eq!(at.state(0), SsrLocusState::Crossed);
        assert_eq!(atg.state(0), SsrLocusState::Crossed);
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
            let mut reads: Vec<u16> = list.iter().map(|d| d.read).collect();
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
        let ssr = ssr_of(&evidence, 0, AT_SIX_REPEATS);

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
        records.walked.set(0);
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

    // ---- what a piece of the census is ---------------------------------------

    /// **The enumeration order is the key's own ordering**, so a fit that walks a sample's
    /// sections is deterministic without anything sorting: ordinary positions before tracts,
    /// read groups by id, strata by their own key. The mutation this catches is the two enum
    /// variants swapping places, which would read every sample's tracts before its ordinary
    /// positions and break the one dependency between the two halves.
    #[test]
    fn sections_enumerate_ordinary_positions_first_then_by_read_group_and_stratum() {
        let period_one = Stratum {
            period: 1,
            reference_repeats: 8,
        };
        let keys = [
            SectionKey::Ssr(ReadGroupId(1), AT_SIX_REPEATS),
            SectionKey::Generic(ReadGroupId(1)),
            SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS),
            SectionKey::Ssr(ReadGroupId(0), period_one),
            SectionKey::Generic(ReadGroupId(0)),
        ];
        let mut sorted = keys;
        sorted.sort_unstable();

        assert_eq!(
            sorted,
            [
                SectionKey::Generic(ReadGroupId(0)),
                SectionKey::Generic(ReadGroupId(1)),
                // Period 1 at 8 repeats sorts before period 2 at 6: the stratum's own order.
                SectionKey::Ssr(ReadGroupId(0), period_one),
                SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS),
                SectionKey::Ssr(ReadGroupId(1), AT_SIX_REPEATS),
            ]
        );
    }

    #[test]
    fn a_section_key_says_whose_reads_and_which_stratum() {
        let generic = SectionKey::Generic(ReadGroupId(3));
        assert_eq!(generic.read_group(), ReadGroupId(3));
        assert_eq!(
            generic.stratum(),
            None,
            "ordinary positions are in no stratum"
        );

        let tracts = SectionKey::Ssr(ReadGroupId(3), AT_SIX_REPEATS);
        assert_eq!(tracts.read_group(), ReadGroupId(3));
        assert_eq!(tracts.stratum(), Some(AT_SIX_REPEATS));
    }

    /// A section filed under a key of the other kind would index one half of the census by the
    /// other's positions, and nothing downstream would notice — so the constructor refuses it.
    #[test]
    #[should_panic(expected = "filed under a key of the other kind")]
    fn a_section_filed_under_the_wrong_kind_of_key_is_refused() {
        let _ = Sections::resident(BTreeMap::from([(
            SectionKey::Generic(ReadGroupId(0)),
            Section::Ssr(SsrEvidence::never_walked(4)),
        )]));
    }

    #[test]
    fn a_resident_sample_lists_its_sections_and_hands_back_the_one_asked_for() {
        let generic = SectionKey::Generic(ReadGroupId(0));
        let tracts = SectionKey::Ssr(ReadGroupId(0), AT_SIX_REPEATS);
        let sections = Sections::resident(BTreeMap::from([
            (generic, Section::Generic(GenericEvidence::never_walked(3))),
            (tracts, Section::Ssr(SsrEvidence::never_walked(2))),
        ]));

        assert_eq!(sections.len(), 2);
        assert_eq!(sections.keys().collect::<Vec<_>>(), vec![generic, tracts]);
        assert!(sections.holds(generic) && sections.holds(tracts));
        assert!(
            !sections.holds(SectionKey::Ssr(ReadGroupId(1), AT_SIX_REPEATS)),
            "a section this sample was never asked for is absent, not empty"
        );
    }

    /// Two sections may never share a byte, since a call reading one seeks to its own offset
    /// and takes its own length.
    #[test]
    fn byte_extents_meet_without_overlapping() {
        let first = ByteExtent::new(0, 100);
        let second = ByteExtent::new(100, 40);
        assert!(!first.overlaps(second), "abutting is not overlapping");
        assert!(first.overlaps(ByteExtent::new(99, 1)));
        assert!(!first.overlaps(ByteExtent::new(100, 0)));
        assert_eq!((second.offset(), second.len()), (100, 40));
        assert!(ByteExtent::new(7, 0).is_empty());
    }

    // ---- the whole cohort ----------------------------------------------------

    /// A resident census never fails to lend, so its calls unwrap in a test.
    fn lent_ok<R>(lent: Result<R, CensusError>) -> R {
        lent.expect("a resident census has no file to fail on")
    }

    /// Two samples over the same one-tract selection, so a cohort can be built from them.
    fn two_sample_cohort() -> Vec<SampleCensusEvidence> {
        ["a", "b"]
            .iter()
            .map(|name| {
                let mut writer = ssr_writer();
                writer.add_locus(&ssr_locus(vec![observation(b"ATATATATATAT", 0, 2)]));
                let mut evidence = writer.finish();
                evidence.sample = (*name).to_string();
                evidence
            })
            .collect()
    }

    /// **The refusal is at the door, before a section is read** (spec §5). Two samples whose
    /// depth ladders differ hold codes meaning different depths while every other value agrees,
    /// and the cohort has to say so by name rather than fit them together.
    #[test]
    fn a_cohort_refuses_two_samples_that_did_not_record_the_same_thing() {
        let mut samples = two_sample_cohort();
        samples[1].terms.depth_ladder = DepthLadderDigest([0; 16]);
        let refusal =
            CohortCensusEvidence::new(samples).expect_err("two ladders is two kinds of evidence");
        assert_eq!(refusal.first, "a");
        assert_eq!(refusal.second, "b");
        assert_eq!(refusal.field, "depth ladder edges");
    }

    /// A cohort lends every sample's section for the band asked for, and **one row a sample**
    /// — the fit reads a stratum from every sample with reads in it at once.
    #[test]
    fn a_cohort_lends_a_band_across_every_sample_and_nothing_outside_it() {
        let mut cohort =
            CohortCensusEvidence::new(two_sample_cohort()).expect("both recorded the same thing");
        assert_eq!(cohort.len(), 2);
        assert_eq!(cohort.read_groups(), &[ReadGroupId(0)]);
        assert_eq!(cohort.strata(), &[AT_SIX_REPEATS]);

        let crossing = lent_ok(cohort.with_strata(&[AT_SIX_REPEATS], |lent| {
            assert_eq!(lent.len(), 2, "one row a sample, not one value");
            lent.iter()
                .map(|sample| {
                    sample
                        .iter()
                        .map(|(_, _, records)| u64::from(records.offsets(0).total()))
                        .sum::<u64>()
                })
                .collect::<Vec<u64>>()
        }));
        assert_eq!(crossing, vec![2, 2], "each sample's two crossing reads");

        // A band naming a stratum this cohort does not hold lends nothing rather than a zero
        // that would read as a sample with no reads there.
        let empty = lent_ok(cohort.with_strata(&[ATG_FOUR_REPEATS], |lent| {
            lent.iter().map(|sample| sample.len()).sum::<usize>()
        }));
        assert_eq!(empty, 0);

        let positions = lent_ok(cohort.with_generic(&[ReadGroupId(0)], |lent| {
            lent.iter().map(|sample| sample.len()).sum::<usize>()
        }));
        assert_eq!(positions, 2, "one ordinary-position section a sample");
    }

    // ---- the twelve values ----------------------------------------------

    /// One read group's ordinary-position section, as a test reads it back from a writer.
    ///
    /// **A test helper and not a public accessor**: what a caller outside this module gets is
    /// the scoped calls, which lend a section and take it back.
    fn generic_of(evidence: &SampleCensusEvidence, group: u32) -> &GenericEvidence {
        evidence
            .generic_sections()
            .into_iter()
            .find(|(id, _)| *id == ReadGroupId(group))
            .map(|(_, records)| records)
            .expect("the read group has an ordinary-position section")
    }

    /// One read group's tracts for one stratum — the same, on the other half.
    fn ssr_of(evidence: &SampleCensusEvidence, group: u32, stratum: Stratum) -> &SsrEvidence {
        evidence
            .ssr_sections()
            .into_iter()
            .find(|(id, at, _)| *id == ReadGroupId(group) && *at == stratum)
            .map(|(_, _, records)| records)
            .expect("the read group has a section for this stratum")
    }

    fn terms() -> RecordingTerms {
        RecordingTerms {
            selection: SelectionTermsDigest::of(&selection_terms()),
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

    /// **The selection travels as a digest a field, and a mismatch is still named** — which is
    /// the whole reason it is a table rather than one digest over the seven.
    #[test]
    fn a_selection_difference_is_named_by_the_selection_rather_than_swallowed() {
        let base = terms();
        let mut moved = selection_terms();
        moved.seed += 1;
        let mut other = base.clone();
        other.selection = SelectionTermsDigest::of(&moved);
        assert_eq!(base.first_disagreement(&other), Some("selection seed"));

        let mut widened = selection_terms();
        widened.generic_target += 1;
        let mut later = base.clone();
        later.selection = SelectionTermsDigest::of(&widened);
        assert_eq!(
            base.first_disagreement(&later),
            Some("generic target position count"),
            "a field further down the list is named as itself, not as the first one"
        );
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
        let (depth, alleles) = generic_of(&evidence, 0).at(0);

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
        let (depth, alleles) = generic_of(&evidence, 0).at(0);

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
        let (depth, alleles) = generic_of(&evidence, 0).at(0);

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

        let generic = generic_of(&records, 0);
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

        let generic = generic_of(&records, 0);
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
        let generic = generic_of(&records, 0);
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

    /// The stratum every tract in these fixtures is in: a two-base motif at six copies.
    const AT_SIX_REPEATS: Stratum = Stratum {
        period: 2,
        reference_repeats: 6,
    };

    /// The second stratum, for the fixtures that need two: a three-base motif at four copies,
    /// so both tracts are twelve bases and only the stratum differs.
    const ATG_FOUR_REPEATS: Stratum = Stratum {
        period: 3,
        reference_repeats: 4,
    };

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

    /// The same, plus a second twelve-base tract at 200 whose motif is three bases — so the two
    /// tracts differ in their stratum and in nothing else.
    fn two_stratum_writer() -> CensusWriter {
        use crate::ng::region_typing::segment_criteria::SsrSegment;
        use crate::ng::repeat_catalog::strata::StratumSampler;

        let at = SsrSegment::new(
            "c0".into(),
            100,
            111,
            crate::ng::types::Motif::new(b"AT").expect("a two-base motif"),
            1.0,
        )
        .expect("a well-formed tract");
        let atg = SsrSegment::new(
            "c0".into(),
            200,
            211,
            crate::ng::types::Motif::new(b"ATG").expect("a three-base motif"),
            1.0,
        )
        .expect("a well-formed tract");
        let mut sampler = StratumSampler::new(16, 0);
        sampler.offer(2, 6, at);
        sampler.offer(3, 4, atg);
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

    /// The three-base-motif tract at 200, which `two_stratum_writer` keeps.
    fn atg_locus(observations: Vec<SequenceObservation>) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(200),
                end: Position(211),
            },
            reference_bases: b"ATGATGATGATG".to_vec().into_boxed_slice(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: LocusKind::Ssr(crate::ng::locus_generation::SsrDetail {
                motif: crate::ng::types::Motif::new(b"ATG").expect("a three-base motif"),
                left_flank: Box::new([]),
                right_flank: Box::new([]),
            }),
        }
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
        let ssr = ssr_of(&records, 0, AT_SIX_REPEATS);

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

    /// **The defect the milestone-A review found, and the one distinction this field exists to
    /// make.** Two reads differing from the reference in *different* places used to come back
    /// as read 0 twice, because the read was numbered within one observation and the walk folds
    /// reads carrying identical bases into one observation. A consumer grouping by read then
    /// saw one read with two interruptions — a bad read — where the truth is two reads with one
    /// each, which is two errors. The fixture beside this one (two reads, *one* interruption)
    /// sits in the single regime where the two numberings agree.
    #[test]
    fn two_reads_carrying_different_interruptions_are_two_reads_and_not_one() {
        let mut writer = ssr_writer();
        writer.add_locus(&ssr_locus(vec![
            observation(b"ATATGTATATAT", 0, 1),
            observation(b"ATATATGTATAT", 0, 1),
        ]));
        let evidence = writer.finish();
        let ssr = ssr_of(&evidence, 0, AT_SIX_REPEATS);

        assert_eq!(
            ssr.differences().len(),
            2,
            "one interruption on each of two reads"
        );
        assert_eq!(
            ssr.differences()[0].offset,
            4,
            "the first read's G sits at the fifth base of the tract"
        );
        assert_eq!(
            ssr.differences()[1].offset,
            6,
            "the second read's at the seventh"
        );
        let mut reads: Vec<u16> = ssr.differences().iter().map(|d| d.read).collect();
        reads.sort_unstable();
        reads.dedup();
        assert_eq!(reads.len(), 2, "two reads, not one read seen twice");
    }

    /// The numbering runs across the whole locus, so a read in the second observation is not
    /// read 0 again — and reads folded into one observation still get one number each.
    #[test]
    fn a_locus_numbers_its_reads_from_zero_once_and_not_once_an_observation() {
        let mut writer = ssr_writer();
        writer.add_locus(&ssr_locus(vec![
            // Two reads sharing an interruption at the fifth base, then one carrying its own
            // at the seventh, then a clean read.
            observation(b"ATATGTATATAT", 0, 2),
            observation(b"ATATATGTATAT", 0, 1),
            observation(b"ATATATATATAT", 0, 4),
        ]));
        let evidence = writer.finish();
        let ssr = ssr_of(&evidence, 0, AT_SIX_REPEATS);

        let at_offset = |offset: i16| {
            let mut reads: Vec<u16> = ssr
                .differences()
                .iter()
                .filter(|d| d.offset == offset)
                .map(|d| d.read)
                .collect();
            reads.sort_unstable();
            reads
        };
        assert_eq!(
            at_offset(4),
            vec![0, 1],
            "the shared interruption, two reads"
        );
        assert_eq!(
            at_offset(6),
            vec![2],
            "the third read of the locus, not the first of its observation"
        );
        assert_eq!(
            ssr.bases_compared(),
            7 * 12,
            "all seven reads were compared, clean ones included"
        );
    }

    /// **The depth end of the range, where one byte failed.** A tract may enter up to 1,000
    /// reads, and while the read number was a byte every read from the 256th onwards was
    /// written as 255 — on the human benchmark trio, at about 265 crossing reads a tract a
    /// sample, that was one mismatch in five. Three hundred reads carrying an interruption
    /// must come back as three hundred reads.
    #[test]
    fn three_hundred_reads_at_one_tract_are_three_hundred_reads() {
        let mut writer = ssr_writer();
        writer.add_locus(&ssr_locus(vec![
            observation(b"ATATGTATATAT", 0, 100),
            observation(b"ATATATGTATAT", 0, 100),
            observation(b"ATATATATGTAT", 0, 100),
        ]));
        let evidence = writer.finish();
        let ssr = ssr_of(&evidence, 0, AT_SIX_REPEATS);

        let mut reads: Vec<u16> = ssr.differences().iter().map(|d| d.read).collect();
        assert_eq!(
            reads.len(),
            300,
            "one interruption on each of three hundred"
        );
        reads.sort_unstable();
        reads.dedup();
        assert_eq!(
            reads,
            (0..300).collect::<Vec<u16>>(),
            "numbered 0 to 299 with none collapsed onto a ceiling"
        );
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
        let ssr = ssr_of(&records, 0, AT_SIX_REPEATS);

        assert_eq!(ssr.offsets(0).at(-1), 3, "one unit short");
        assert_eq!(ssr.offsets(0).at(0), 5);
        assert_eq!(
            ssr.bases_compared(),
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
            generic_of(&records, 0).at(0).0,
            DepthCode::Binned(edges.bin_for(6))
        );
        assert_eq!(
            generic_of(&records, 1).at(0).0,
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
            generic_of(&records, 0).at(0).0,
            DepthCode::Binned(edges.bin_for(5))
        );
        assert_eq!(
            generic_of(&records, 1).at(0).0,
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

        let generic = generic_of(&records, 0);
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
        assert_eq!(generic_of(&records, 0).at(0).0, DepthCode::NeverWalked);
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
