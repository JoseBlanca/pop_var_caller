//! Shared ng vocabulary — the domain newtypes cross-step code speaks. It starts as this
//! one file and splits into concept modules (`units`, `locus`, …) as clusters grow
//! (`doc/devel/ng/arch/module_layout.md` principle 3). Seeded here with only what the
//! `RefSeq` reference accessor needs.

use std::fmt;

/// Which reference sequence a coordinate refers to: an index into the reference contig
/// table ([`crate::fasta::ContigList`]), in `@SQ` / `.fai` order. Unconstrained — any
/// `u32` is a legal index at the type level, and an out-of-range id is caught at fetch
/// time — so the field is public and there is no checked constructor.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct ContigId(pub u32);

impl ContigId {
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

/// A **1-based** reference position — the coordinate ng speaks everywhere
/// (`ng_step_interfaces.md` §1), matching VCF/SAM/IGV and the production engine.
/// Unconstrained: any `u64` is a legal value at the type level, so the field is
/// public and there is no checked constructor.
///
/// The base is the point. ng chose 1-based so that a coordinate printed in a log,
/// a VCF, or a bug report means the same thing everywhere — no mental `+ 1`, and
/// no class of off-by-one that only shows up as a wrong genotype. The exceptions
/// are named and local: `RepeatInterval` is 0-based because it indexes a byte
/// slice, and `regions::RegionSet` is production's (converted at
/// `GenomeRegions`).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Position(pub u64);

impl Position {
    #[inline]
    pub fn get(self) -> u64 {
        self.0
    }
}

/// One base, genome-wide: which contig, and where along it. A [`Position`] on
/// its own does not identify a base — position 1000 exists on every contig —
/// which is why so many signatures thread `(contig, position)` as two
/// parameters. This is that pair given a name, and the point-shaped sibling of
/// [`GenomeRegion`].
///
/// **The derived [`Ord`] is genome order**: contig index first, then position
/// within it, because that is the field order. That is what lets the type serve
/// directly as a sort key wherever reads or loci are ordered along the
/// reference — the read-order guard and the k-way merge are its first uses
/// (`doc/devel/ng/arch/sample_reads.md` §1.1).
///
/// The ordering is only meaningful because every alignment file's `ref_id` was
/// proved equal to its [`ContigId`] when the file was opened
/// (`doc/devel/ng/spec/alignment_file.md` §3.1). Without that gate, contig
/// indices from different files would not be comparable at all.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct GenomePosition {
    pub contig: ContigId,
    pub position: Position,
}

/// A **physical** piece of DNA: a contig plus a 1-based **inclusive** range
/// `[start, end]`. No genetic claim — that is what distinguishes it from a
/// *locus* (`typed_regions.md` §1.1).
///
/// The consolidation `ng_step_interfaces.md` §1 reserved, and the typed-region
/// generator is its first real use: it replaces `regions::Region` (0-based,
/// `u32`, production's) and `bam::ContigInterval` for ng's purposes.
///
/// **Inclusive, not half-open**, and deliberately: a half-open `end` cannot name
/// the last base of a contig without arithmetic, and "the region ends at 21"
/// meaning base 20 is exactly the ambiguity the 1-based decision exists to
/// delete. The cost is one `+ 1` in [`Self::len`], which is where it belongs —
/// stated once, in the type.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct GenomeRegion {
    pub contig: ContigId,
    pub start: Position,
    pub end: Position,
}

impl GenomeRegion {
    /// The region's length in bases. **The one place the inclusive `+ 1` lives.**
    ///
    /// Saturating rather than panicking on an inverted region: this is a plain
    /// data type with public fields (no constructor to enforce `start <= end`),
    /// and a length of 0 is a truer answer for an empty span than a panic in a
    /// getter would be. Callers that require well-formedness say so themselves.
    #[inline]
    pub fn len(self) -> u64 {
        (self.end.get() + 1).saturating_sub(self.start.get())
    }

    /// Whether the region covers no bases (`end < start`).
    #[inline]
    pub fn is_empty(self) -> bool {
        self.end.get() < self.start.get()
    }

    /// Whether `pos` falls inside the region, bounds included.
    #[inline]
    pub fn contains(self, pos: Position) -> bool {
        self.start <= pos && pos <= self.end
    }
}

/// `contig 3:940-1100` — the shape an error message wants.
///
/// **It says `contig 3` and not `chr4` because the type holds an id, not a
/// name.** Rendering the id as though it were a name is a lossy translation this
/// codebase has made before and recorded as a defect
/// (`locus_generation_pileup.md` — A0 deleted one), so the word `contig` is
/// there to stop a reader taking the number for a chromosome. A caller holding
/// the reference's contig table ([`ContigTable`](crate::ng::ref_seq::ContigTable))
/// can do better and should.
impl std::fmt::Display for GenomeRegion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "contig {}:{}-{}",
            self.contig.get(),
            self.start.get(),
            self.end.get()
        )
    }
}

// ---------------------------------------------------------------------
// Scalar newtypes — the domain quantities cross-step code speaks. Seeded
// here with only the scalars read filtering (ng step 1) touches; the rest
// of the `ng_step_interfaces.md` §1 vocabulary lands as later steps need
// it. See `doc/devel/ng/arch/read_filtering.md` §2.1.
// ---------------------------------------------------------------------

/// SAM mapping quality (MAPQ): the aligner's Phred-scaled confidence that the
/// read is placed at the right locus. `0` = "could be anywhere"; `60` = "as
/// sure as this aligner gets". MAPQ unavailable (SAM `0xFF`) is treated as `0`
/// by callers. Unconstrained — every `u8` is a legal value, so the field is
/// public and there is no checked constructor.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct MapQual(pub u8);

impl MapQual {
    #[inline]
    pub fn get(self) -> u8 {
        self.0
    }
}

/// A single base call's Phred quality (0–93). Unconstrained — any `u8` is a
/// legal value, so the field is public and there is no checked constructor.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct BaseQual(pub u8);

impl BaseQual {
    #[inline]
    pub fn get(self) -> u8 {
        self.0
    }
}

/// A length in base pairs — the generic length currency both the SNP/indel and
/// STR paths speak (only *repeat-unit* quantities carry the `Ssr` prefix). Here
/// it measures a read's decoded length. Unconstrained — any `u64` is a legal
/// value, so the field is public and there is no checked constructor.
///
/// `u64` since B2 (spec §4): ng speaks one width, so nothing narrows, nothing is
/// checked, and no off-by-width bug is possible. Ids stay `u32` — they index a
/// table, they are not positions.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Bp(pub u64);

impl Bp {
    #[inline]
    pub fn get(self) -> u64 {
        self.0
    }
}

/// Which read group a read came from: an index into the run's read-group table
/// (`ReadGroups`, `src/ng/read/input/read_groups.rs`).
///
/// A **read group** is one SAM `@RG` record — in practice one library preparation
/// sequenced in one run. It is the unit a per-chemistry error model keys on,
/// because PCR stutter and per-base error are properties of the library
/// preparation and of the DNA's condition, not of the individual the DNA came
/// from (`doc/devel/ng/spec/read_groups.md` §1).
///
/// Ids are minted only when the table is built, in input-file order and then
/// header order within a file, so the same input list always yields the same ids
/// however the files were read (spec §4). Unconstrained — any `u32` is a legal
/// index at the type level, and an out-of-range id is caught at lookup — so the
/// field is public and there is no checked constructor. `u32` rather than `u64`
/// for the reason [`Bp`] gives: an id indexes a table, it is not a position.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct ReadGroupId(pub u32);

impl ReadGroupId {
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

/// Just the index — `3`, not `ReadGroupId(3)`. A message naming a read group supplies its
/// own word for it ("read group {read_group}"), so the type renders the number alone.
///
/// **It renders an index into the run's read-group table, not a library name.** A caller
/// holding that table can print the `@RG ID` and should; this is what an error message can
/// say without one.
impl fmt::Display for ReadGroupId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A probability held as its natural logarithm — the number stored is `ln(p)`, not
/// `p` itself.
///
/// `f64::NEG_INFINITY` is a legal value, not an error: it is `ln(0)`, the score of
/// something impossible — a read line-up that cannot happen. That is the whole reason
/// probabilities are carried this way here: an impossible line-up reaches a finite
/// sentinel a caller can see (`-∞`), where in ordinary probabilities it reaches `0`,
/// indistinguishable from a value that merely got too small to represent. Every finite
/// `f64` and `-∞` is a valid log-probability, so — like [`Bp`], and unlike
/// [`MismatchFraction`] — the value is unconstrained, the field is public, and there is
/// no checked constructor.
///
/// Its point as a *distinct type* is that the compiler refuses to mix it with an
/// ordinary (linear) probability. That is the mistake the alignment module is most
/// exposed to: the production code it ports returns linear probabilities, the
/// conversion to a logarithm happens at one boundary, and a `LogProb` accidentally
/// handed a raw probability would be a plausible wrong number rather than a compile
/// error without this wrapper.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct LogProb(pub f64);

impl LogProb {
    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// A fraction of mismatched bases, constrained to `[0, 1]`. Unlike the
/// unconstrained newtypes above, an out-of-range value is *unrepresentable*:
/// the field is private and construction goes through the checked
/// [`Self::try_new`]. Read filtering uses it as the mismatch-fraction threshold
/// (filter #8), whose source is an untrusted CLI/config value — so the policy
/// is fail loudly, never silently coerce.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct MismatchFraction(f32);

impl MismatchFraction {
    /// The only constructor. A fraction outside `[0, 1]` is a user error —
    /// reject it rather than coerce.
    pub fn try_new(x: f32) -> Result<Self, DomainError> {
        (0.0..=1.0)
            .contains(&x)
            .then_some(Self(x))
            .ok_or(DomainError::MismatchFraction(x))
    }

    #[inline]
    pub fn get(self) -> f32 {
        self.0
    }
}

// ---------------------------------------------------------------------
// The parameters a caller runs on — three scalars step 4 fits (the error
// rate, the genotype frequencies, the inbreeding coefficient) and one it
// is handed (ploidy). Four types and not one shared `Probability`: three
// of them are fractions in `[0, 1]`, so a single type would let an
// inbreeding coefficient be handed to something expecting an error rate
// and compile (`arch/parameter_prepass_generic.md` §2.1).
//
// They live in the shared vocabulary rather than in `parameter_estimation/`
// because their consumers are *other* steps: the likelihood (step 7) reads
// the error rate, the genotype prior (step 8) reads the genotype
// frequencies and the inbreeding coefficient, and ploidy reaches both.
// **Neither step exists yet** — defining the types here is what keeps them
// from importing out of a sibling step's module when they arrive. They are
// also the seed of the `genotype`/`params` split `module_layout.md`
// principle 3 anticipates for this file.
//
// Each follows `MismatchFraction`'s shape: private field, checked
// `try_new`, `.get()`. `try_new` is the **boundary** constructor — it
// rejects rather than coerces, for values arriving from outside the
// program. The fits construct through the same door and `.expect()`,
// because a frequency off the simplex means our own arithmetic is broken
// and there is nothing a caller could do about it.
// ---------------------------------------------------------------------

/// `Ok(x)` when `x` is a probability — inside `[0, 1]`, both endpoints included.
/// Otherwise the caller's own [`DomainError`] variant, passed in as its tuple
/// constructor so the message still names the quantity that was wrong.
///
/// **One predicate, written once.** `NaN`, `+∞` and `-∞` all fail it without any
/// help from `is_finite`: `contains` is `0.0 <= x && x <= 1.0`, no comparison with
/// `NaN` is true, `+∞` is not `<= 1` and `-∞` is not `>= 0`. Constructors
/// spelling the same range test separately is how one of them ends up written
/// `0.0..1.0` and rejecting a genotype frequency of exactly one — a real answer for
/// a fully homozygous sample.
///
/// `pub(crate)` for that same reason rather than for convenience: the STR path's three
/// slippage rates are constrained the same way
/// (`parameter_estimation::ssr::slippage`), so **six** constructors now share this
/// predicate. It is not hypothetical drift — `SiteNoise::try_new`
/// (`parameter_estimation::generic`) already spells the range test by hand, which is the
/// seventh probability in the crate and the one this predicate did not reach.
pub(crate) fn checked_probability(
    x: f64,
    reject: fn(f64) -> DomainError,
) -> Result<f64, DomainError> {
    if (0.0..=1.0).contains(&x) {
        Ok(x)
    } else {
        Err(reject(x))
    }
}

/// A per-base sequencing error rate: how often a read shows a base other than the
/// one on the template it was read from. A probability in `[0, 1]`.
///
/// Estimated **per read group**, because the chemistry belongs to the library
/// preparation and the sequencing run, not to the individual whose DNA they read
/// (`spec/parameter_prepass_generic.md` §2).
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct ErrorRate(f64);

impl ErrorRate {
    /// The only constructor. A rate that is not a probability in `[0, 1]` is
    /// rejected rather than coerced.
    pub fn try_new(rate: f64) -> Result<Self, DomainError> {
        checked_probability(rate, DomainError::ErrorRate).map(Self)
    }

    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// How common one genotype is in a sample's genome — the share of sites carrying
/// it. A probability in `[0, 1]`; a set of them over the genotypes at one ploidy
/// sums to one.
///
/// On the diploid path the three are the homozygous-reference rate, the
/// heterozygosity, and the homozygous-non-reference rate.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct GenotypeFrequency(f64);

impl GenotypeFrequency {
    /// The only constructor. A frequency that is not a probability in `[0, 1]` is
    /// rejected rather than coerced.
    pub fn try_new(frequency: f64) -> Result<Self, DomainError> {
        checked_probability(frequency, DomainError::GenotypeFrequency).map(Self)
    }

    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// The inbreeding coefficient: the fraction of an individual's analysable genome
/// lying in runs of homozygosity, where the two copies descend from one recent
/// ancestor. `0` for an outcrosser, approaching `1` for a long-selfed line.
///
/// A user may supply one on the command line, which is why the constructor
/// rejects rather than coerces (`spec/parameter_prepass_generic.md` §6.4).
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct InbreedingF(f64);

impl InbreedingF {
    /// The only constructor. A coefficient that is not a fraction in `[0, 1]` is
    /// rejected rather than coerced.
    pub fn try_new(coefficient: f64) -> Result<Self, DomainError> {
        checked_probability(coefficient, DomainError::InbreedingF).map(Self)
    }

    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// How many copies of the genome an individual carries at a region — two on a
/// diploid autosome, one on a haploid sex chromosome.
///
/// An input to every fit rather than a global constant: it varies by region
/// within one genome (`spec/parameter_prepass.md` §3).
///
/// **Constrained, unlike the unconstrained newtypes elsewhere in this file**
/// ([`ContigId`], [`Position`], [`Bp`], [`LogProb`]): the likelihood divides by the
/// number of copies, so a zero is a division by zero rather than an odd answer.
/// `Ord` because it keys the histogram and output maps, where the derived order is
/// the natural one — fewest copies first.
///
/// No upper bound. A polyploid crop is in scope — ploidy varies by region within
/// one genome — so any ceiling short of `u8::MAX` would reject a real hexaploid or
/// octoploid region.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Ploidy(u8);

impl Ploidy {
    /// The only constructor. Zero copies is rejected: it is not a genome.
    pub fn try_new(copies: u8) -> Result<Self, DomainError> {
        if copies == 0 {
            return Err(DomainError::Ploidy(copies));
        }
        Ok(Self(copies))
    }

    #[inline]
    pub fn get(self) -> u8 {
        self.0
    }
}

/// Just the copy number — `2`, not `Ploidy(2)`. An error message that names a ploidy
/// supplies its own word for it ("at ploidy {ploidy}"), so the type renders the number
/// and nothing else.
impl fmt::Display for Ploidy {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// A domain-invariant violation — the ng-wide error raised when an untrusted
/// value falls outside a constrained newtype's range. Introduced with its
/// first variant; later constrained types (`AlleleFreq`, `Theta`, …) add their
/// own variants as they arrive. `#[non_exhaustive]` so matchers accept those
/// future variants without breaking.
///
/// **`PartialEq` is IEEE equality on the `f64` payloads, so an error carrying a
/// `NaN` is not equal to itself.** A `NaN` input is exactly what the three rate
/// constructors reject, so this is not a corner case: compare such a rejection with
/// `matches!(err, Err(DomainError::ErrorRate(r)) if r.is_nan())`, never with
/// `assert_eq!`, which fails printing two sides that render identically.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum DomainError {
    /// A [`MismatchFraction`] was constructed from a value outside `[0, 1]`.
    #[error("mismatch fraction {0} is outside [0, 1]")]
    MismatchFraction(f32),

    /// A `SiteNoise` was constructed with a noisy-site share outside `[0, 1]`.
    ///
    /// Its own variant rather than a shared "not a probability", for the reason the four
    /// scalars above have four: a share of *sites* and a rate per *read* are different
    /// quantities, and a message naming the wrong one sends the reader to the wrong fit.
    #[error("noisy-site fraction {0} is outside [0, 1]")]
    SiteNoiseFraction(f64),
    /// A per-base error rate is not a probability in `[0, 1]`.
    ///
    /// **Three constructors raise this one variant** — [`ErrorRate::try_new`],
    /// `FlatEmission::try_new` and `SsrSequenceMarginal::try_new` — so the message
    /// names the quantity but not which of them rejected the value. Deliberate
    /// reuse rather than an oversight: all three mean the same thing. If a log line
    /// ever needs to tell them apart, that is a variant split, not a message edit.
    #[error("per-base error rate {0} is not a finite probability in [0, 1]")]
    ErrorRate(f64),
    /// A [`GenotypeFrequency`] was built from a value that is not a finite
    /// probability in `[0, 1]`.
    #[error("genotype frequency {0} is not a finite probability in [0, 1]")]
    GenotypeFrequency(f64),
    /// An [`InbreedingF`] was built from a value that is not a finite fraction in
    /// `[0, 1]`.
    #[error("inbreeding coefficient {0} is not a finite fraction in [0, 1]")]
    InbreedingF(f64),
    /// A [`Ploidy`] was built from zero genome copies, which the likelihood
    /// divides by.
    #[error("ploidy {0} is not a positive number of genome copies")]
    Ploidy(u8),
    /// A read's bases and its qualities were paired but differ in length.
    #[error("read has {bases} bases but {qualities} qualities")]
    ReadQualityLengthMismatch { bases: usize, qualities: usize },

    /// An [`SsrPeriod`] was built from zero bases, or from more than a repeat
    /// unit can hold.
    ///
    /// **Zero is the one that has to be unrepresentable.** A tract's length
    /// becomes a repeat count by dividing by its period, so a period of zero is
    /// a division by zero rather than an odd answer
    /// (`arch/parameter_prepass_ssr.md` §2.1).
    ///
    /// Carries a `usize` rather than the `u8` the period is stored in, so that
    /// the message names the value the **caller** offered: every producer of a
    /// period in this crate holds a `usize` ([`Motif::period`], a tract length
    /// divided by one), and a variant narrower than its callers is what makes
    /// `try_new(n as u8)` the natural call — under which 258 arrives as 2 and
    /// validates as a dinucleotide.
    #[error("repeat period {0} is outside the STR period range 1..={MAX_MOTIF_LEN}")]
    SsrPeriod(usize),

    /// A slippage rate — how often a read at a repeat tract shows a length
    /// other than its allele's — is not a probability in `[0, 1]`.
    #[error("slippage rate {0} is not a finite probability in [0, 1]")]
    SlipRate(f64),

    /// The share of slipped reads that **gained** repeats rather than losing
    /// them is not a probability in `[0, 1]`.
    ///
    /// Its own variant beside [`Self::SlipRate`] and [`Self::SlipStepDecay`],
    /// for the reason the four scalars above have four: all three are fractions
    /// in `[0, 1]`, and a message naming the wrong one sends a reader to the
    /// wrong parameter of the same fit.
    ///
    /// **The message says what the share is of**, because the two readings
    /// differ by a factor of the slippage rate — about fiftyfold at a level of
    /// 0.02 — and a bare "gain share" reads either way.
    #[error("gain share of slipped reads {0} is not a finite probability in [0, 1]")]
    SlipGainShare(f64),

    /// The chance that a slipped read moved a second repeat, given that it
    /// moved a first, is not a probability in `[0, 1]`.
    #[error("slip step decay {0} is not a finite probability in [0, 1]")]
    SlipStepDecay(f64),
}

// ---------------------------------------------------------------------
// The motif — STR domain vocabulary, shared across steps
// ---------------------------------------------------------------------

/// STR scope: a repeat unit (period) is between 1 and this many bases.
pub const MAX_MOTIF_LEN: usize = 6;

/// A motif's bytes were not a valid STR period.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[non_exhaustive]
pub enum MotifError {
    /// Length is `0` or above [`MAX_MOTIF_LEN`] — outside the STR period range.
    #[error("motif length {len} is outside the STR period range 1..={MAX_MOTIF_LEN}")]
    BadLength { len: usize },
}

/// A tandem-repeat unit — the repeat's period, 1..=[`MAX_MOTIF_LEN`] bases.
///
/// The bytes are stored **verbatim**: the reference-strand, phase-faithful unit
/// exactly as it tiles the locus (e.g. `CAG`), *not* canonicalized. Rotating to
/// a canonical form (`CAG` → `AGC`) would break tiling, and reconstruction reads
/// phase-correct bytes off the reference anyway; the canonical *class* used for
/// stutter pooling is derived on demand downstream, never stored here.
///
/// Inline and `Copy`: a fixed 6-byte buffer plus a length, never heap-allocated.
/// Unused tail bytes are zero, so the derived `Eq`/`Hash` compare only the live
/// prefix (`0` is not a valid base, so it cannot collide with one).
///
/// INVARIANT (relied on by the derived `Eq`/`Hash`): every constructor MUST
/// zero-initialize the unused tail of `buf`. [`Motif::new`] is currently the
/// only one; any future constructor must uphold this or the derived impls will
/// treat equal motifs as distinct.
///
/// ## Why this is a port, when spec §4 said to reuse production's
///
/// Spec §4 kept `ssr::types::Motif` on the grounds that it *"carries no
/// coordinates and no width, and so has nothing to rebase — the Revision's
/// 'reuse where it costs production nothing' case exactly."* The first half is
/// true; **the conclusion was wrong**, and the compiler said so.
///
/// `ssr::types::Motif` is `pub(crate)`. ng's [`SsrSegment`] is `pub` (the ng-sibling
/// convention) and returns a motif, so reusing it trips rustc's
/// `private_interfaces` lint — a `pub` item leaking a `pub(crate)` type. The
/// three ways out: widen `Motif` in `src/ssr/types.rs` (**touching production —
/// forbidden**); demote ng's whole classification surface to `pub(crate)` (bends ng's
/// convention *and* buys `dead_code` warnings for every item until its Milestone
/// D consumer exists); or port the 40 lines. So reuse did **not** cost
/// production nothing — it cost a visibility compromise, which is precisely the
/// coupling "a fresh ng caller from scratch" (owner, 2026-07-16) exists to
/// avoid.
///
/// Ported, ng's `region_typing` names nothing from `src/ssr/` outside its
/// `#[cfg(test)]` differential — which is exactly where the dependency belongs.
/// The type is coordinate-free and trivially checkable, so the duplication is
/// cheap and the drift risk is near zero; the differential compares motifs by
/// bytes and would catch it anyway.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
pub struct Motif {
    buf: [u8; MAX_MOTIF_LEN],
    len: u8,
}

impl Motif {
    /// Build a motif from its bytes, validating the STR period range.
    ///
    /// The bytes are taken verbatim (no canonicalization); the caller is
    /// responsible for supplying the phase-faithful, reference-strand unit.
    pub fn new(bytes: &[u8]) -> Result<Self, MotifError> {
        let len = bytes.len();
        if len == 0 || len > MAX_MOTIF_LEN {
            return Err(MotifError::BadLength { len });
        }
        let mut buf = [0u8; MAX_MOTIF_LEN];
        buf[..len].copy_from_slice(bytes);
        Ok(Self {
            buf,
            len: len as u8,
        })
    }

    /// The motif bytes.
    #[inline]
    pub fn as_bytes(&self) -> &[u8] {
        &self.buf[..self.len as usize]
    }

    /// The period (motif length, in bases).
    #[inline]
    pub fn period(&self) -> usize {
        self.len as usize
    }

    /// The period as the checked domain type, for a consumer that divides by it.
    ///
    /// **Beside [`Self::period`] rather than replacing it.** That accessor returns a
    /// `usize` and has callers across three modules; changing its type under them to
    /// serve one new consumer is a change to code that is working. A motif is
    /// constructed only through [`Self::new`], which already rejects a length outside
    /// `1..=MAX_MOTIF_LEN`, so this conversion cannot fail — the `expect` is a claim
    /// about that invariant, not a fallible path.
    #[inline]
    pub fn ssr_period(&self) -> SsrPeriod {
        SsrPeriod::try_new(self.period()).expect("a motif's length is checked at construction")
    }
}

/// A repeat unit's length in bases — the **period** of a tandem repeat, `1..=MAX_MOTIF_LEN`.
///
/// **The first half of the axis a stutter model is stratified by** (`spec/parameter_prepass_ssr.md`
/// §4): how much a tract slips depends on its motif's period and on how many copies of that motif
/// it holds, and those two together name the group of loci one set of stutter parameters is fitted
/// from.
///
/// **Constrained, unlike [`Motif::period`]'s `usize`**: a tract's length becomes a repeat count by
/// dividing by the period, so zero is a division by zero. The upper bound is the STR scope every
/// other part of ng already works in ([`MAX_MOTIF_LEN`]).
///
/// Here rather than in the step that fits stutter, because the likelihood (step 7) and the genotype
/// prior (step 8) both read a locus's stutter parameters and so both name the period they were
/// fitted at (`arch/module_layout.md`).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct SsrPeriod(u8);

impl SsrPeriod {
    /// The only constructor. A period outside `1..=MAX_MOTIF_LEN` is rejected rather
    /// than coerced.
    ///
    /// **Takes a `usize`, though it stores a `u8`, and the width is the point.** Every
    /// producer of a period in this crate holds a `usize` — [`Motif::period`], and a
    /// tract's length divided by one — so a `u8` parameter would make `try_new(n as u8)`
    /// the natural call at each of them, and that cast turns 258 into 2: a rejected value
    /// arriving as an accepted dinucleotide. Widening the door is what keeps the check
    /// meaningful at the call sites that actually exist.
    pub fn try_new(bases: usize) -> Result<Self, DomainError> {
        if bases == 0 || bases > MAX_MOTIF_LEN {
            return Err(DomainError::SsrPeriod(bases));
        }
        // Inside `1..=MAX_MOTIF_LEN`, which is 6, so the narrowing cannot lose anything.
        Ok(Self(bases as u8))
    }

    #[inline]
    pub fn get(self) -> u8 {
        self.0
    }
}

/// Just the number of bases — `3`, not `SsrPeriod(3)`. A message naming a period supplies
/// its own word for it ("period {period}"), so the type renders the number and nothing else.
impl fmt::Display for SsrPeriod {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl fmt::Debug for Motif {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Bases are ASCII; render as text for readable test output / logs,
        // falling back to bytes if a motif ever held non-UTF-8.
        match std::str::from_utf8(self.as_bytes()) {
            Ok(s) => write!(f, "Motif({s:?})"),
            Err(_) => write!(f, "Motif({:?})", self.as_bytes()),
        }
    }
}

// Where `Motif` used to live, and why it moved. It was written in
// `region_typing::segment_criteria`, because that is where the STR classifier needed it. It
// then gained consumers in `locus_generation` and `alignment` — three modules across
// pipeline-stage boundaries — and `alignment` states in its own module doc that it is **not
// a pipeline step and knows no callers**, which an import from step 3 flatly contradicted,
// on ng's *public* surface. `module_layout.md` principle 3 already assigns shared STR domain
// vocabulary here. `segment_criteria` re-exports it, so step 3's own callers are untouched.

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mismatch_fraction_accepts_boundary_values() {
        assert_eq!(MismatchFraction::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(MismatchFraction::try_new(1.0).unwrap().get(), 1.0);
        assert_eq!(MismatchFraction::try_new(0.10).unwrap().get(), 0.10);
    }

    #[test]
    fn mismatch_fraction_rejects_out_of_range() {
        assert_eq!(
            MismatchFraction::try_new(-0.01),
            Err(DomainError::MismatchFraction(-0.01))
        );
        assert!(matches!(
            MismatchFraction::try_new(1.01),
            Err(DomainError::MismatchFraction(_))
        ));
        // The `[0, 1]` range check rejects the infinities as well.
        assert!(MismatchFraction::try_new(f32::INFINITY).is_err());
        assert!(MismatchFraction::try_new(f32::NEG_INFINITY).is_err());
    }

    #[test]
    fn mismatch_fraction_rejects_nan() {
        // NaN is neither <= 1 nor >= 0, so the range check rejects it.
        assert!(MismatchFraction::try_new(f32::NAN).is_err());
    }

    #[test]
    fn unconstrained_newtypes_expose_their_value() {
        assert_eq!(MapQual(20).get(), 20);
        assert_eq!(BaseQual(93).get(), 93);
        assert_eq!(Bp(150).get(), 150);
        assert_eq!(Position(1).get(), 1);
        assert_eq!(ReadGroupId(7).get(), 7);
    }

    /// `LogProb` is unconstrained: any finite logarithm round-trips, and `-∞` — the
    /// score of an impossible line-up — is a value it must carry, not reject. That
    /// `-∞` is preserved rather than coerced is the property callers rely on to tell
    /// "cannot happen" from "too small to represent".
    #[test]
    fn log_prob_carries_any_logarithm_including_negative_infinity() {
        assert_eq!(LogProb(0.0).get(), 0.0); // ln(1): certainty
        assert_eq!(LogProb(-2.5).get(), -2.5);
        assert_eq!(LogProb(f64::NEG_INFINITY).get(), f64::NEG_INFINITY);
        // Ordering places the impossible line-up below every finite score, which is
        // what lets a caller compare log-probabilities directly.
        assert!(LogProb(f64::NEG_INFINITY) < LogProb(-1000.0));
        assert!(LogProb(-2.5) < LogProb(0.0));
    }

    /// **Why the type derives `PartialOrd`, not `Ord`.** A `NaN` is outside the
    /// documented domain (finite ∪ {−∞}), but the unconstrained public field makes it
    /// *representable*, and a `NaN` compares unordered to everything — including itself.
    /// A total order would be a lie here, and this pins that: a stray `NaN` reaching a
    /// caller that maxes log-scores silently loses every comparison rather than
    /// panicking, so this documents the hazard at the type rather than leaving it to be
    /// rediscovered.
    #[test]
    fn log_prob_partialord_is_not_a_total_order_for_nan() {
        let nan = LogProb(f64::NAN);
        let also_nan = LogProb(f64::NAN);
        let finite = LogProb(0.0);
        // `None` for two NaNs means unordered *and* unequal — reflexivity fails, so
        // `Eq`/`Ord` would be unsound and are correctly not derived.
        assert_eq!(nan.partial_cmp(&also_nan), None);
        // A NaN is unordered against every finite score, in either direction — the
        // silent-comparison hazard a caller that maxes log-scores must know about.
        assert_eq!(nan.partial_cmp(&finite), None);
        assert_eq!(finite.partial_cmp(&nan), None);
        assert!(nan != finite);
    }

    /// The field is unconstrained, so `+∞` — an out-of-domain value the doc does not
    /// name as valid — is carried verbatim rather than coerced or rejected, exactly as
    /// `-∞` is. The type is a transparent wrapper: what goes in comes back out.
    #[test]
    fn log_prob_carries_positive_infinity_out_of_domain() {
        assert_eq!(LogProb(f64::INFINITY).get(), f64::INFINITY);
        assert!(LogProb(f64::INFINITY) > LogProb(1e300));
    }

    /// The three `[0, 1]` rates accept both endpoints — a genotype frequency of
    /// exactly zero and an inbreeding coefficient of exactly one are both real
    /// answers, so a half-open check would reject valid data.
    #[test]
    fn the_constrained_rates_accept_both_endpoints() {
        assert_eq!(ErrorRate::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(ErrorRate::try_new(1.0).unwrap().get(), 1.0);
        assert_eq!(ErrorRate::try_new(0.001).unwrap().get(), 0.001);

        assert_eq!(GenotypeFrequency::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(GenotypeFrequency::try_new(1.0).unwrap().get(), 1.0);

        assert_eq!(InbreedingF::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(InbreedingF::try_new(1.0).unwrap().get(), 1.0);
    }

    /// Each rate names **its own quantity** when it rejects, so a message cannot
    /// send a reader to the wrong fit. `GenotypeFrequency` and `InbreedingF` have
    /// their own variants; `ErrorRate` shares `DomainError::ErrorRate` with the two
    /// emission models, which is deliberate reuse — all three mean "a per-base
    /// error rate that is not a probability".
    ///
    /// **Both directions, for all three.** Each range check is two comparisons, and
    /// a test that only ever crosses one of them leaves the other free to be
    /// widened: `InbreedingF` accepting `1.5` is the live hazard, since a user
    /// types that one at a shell.
    #[test]
    fn each_constrained_rate_rejects_out_of_range_in_both_directions() {
        assert_eq!(
            ErrorRate::try_new(-0.01),
            Err(DomainError::ErrorRate(-0.01))
        );
        assert_eq!(ErrorRate::try_new(1.01), Err(DomainError::ErrorRate(1.01)));
        assert_eq!(
            GenotypeFrequency::try_new(-0.5),
            Err(DomainError::GenotypeFrequency(-0.5))
        );
        assert_eq!(
            GenotypeFrequency::try_new(1.5),
            Err(DomainError::GenotypeFrequency(1.5))
        );
        assert_eq!(
            InbreedingF::try_new(-0.5),
            Err(DomainError::InbreedingF(-0.5))
        );
        assert_eq!(
            InbreedingF::try_new(1.5),
            Err(DomainError::InbreedingF(1.5))
        );
    }

    /// `NaN` and both infinities are not probabilities and none of them
    /// constructs — so a rate arriving from a division by zero cannot enter a
    /// likelihood.
    ///
    /// **The `[0, 1]` range check rejects all three on its own.** `contains` is
    /// `0.0 <= x && x <= 1.0`: no comparison with `NaN` is true, `+∞` is not `<= 1`,
    /// and `-∞` is not `>= 0`. There is no `is_finite` call anywhere in
    /// `checked_probability`, and none is needed.
    ///
    /// `is_err()` rather than `assert_eq!` on purpose: `DomainError` compares its
    /// `f64` payloads by IEEE equality, under which a `NaN` rejection is not equal
    /// to itself.
    #[test]
    fn the_constrained_rates_reject_nan_and_the_infinities() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(ErrorRate::try_new(bad).is_err(), "error rate {bad}");
            assert!(
                GenotypeFrequency::try_new(bad).is_err(),
                "genotype frequency {bad}"
            );
            assert!(InbreedingF::try_new(bad).is_err(), "inbreeding {bad}");
        }
    }

    /// Ploidy zero is the one value the type exists to make unrepresentable: the
    /// likelihood divides by the number of copies.
    #[test]
    fn ploidy_rejects_zero_and_accepts_every_real_copy_number() {
        assert_eq!(Ploidy::try_new(0), Err(DomainError::Ploidy(0)));
        assert_eq!(Ploidy::try_new(1).unwrap().get(), 1);
        assert_eq!(Ploidy::try_new(2).unwrap().get(), 2);
        assert_eq!(Ploidy::try_new(4).unwrap().get(), 4);
    }

    /// `Ploidy` keys the histogram and the emitted rate maps, so its order has to
    /// be the natural one — a haploid region's cells sort before a diploid's.
    /// `Ploidy` renders as the bare copy number, because the messages that name one
    /// supply their own word for it ("at ploidy {ploidy}"). Tested because nothing else
    /// reads the impl: the error message that uses it asserts the sample and the site
    /// count, so replacing this body with a constant would leave the suite green.
    #[test]
    fn ploidy_displays_as_the_bare_copy_number() {
        assert_eq!(Ploidy::try_new(2).unwrap().to_string(), "2");
        assert_eq!(Ploidy::try_new(1).unwrap().to_string(), "1");
        assert_eq!(Ploidy::try_new(u8::MAX).unwrap().to_string(), "255");
        assert_eq!(
            format!("at ploidy {}", Ploidy::try_new(4).unwrap()),
            "at ploidy 4"
        );
    }

    #[test]
    fn ploidy_orders_by_copy_number() {
        let mut ploidies = [
            Ploidy::try_new(4).unwrap(),
            Ploidy::try_new(1).unwrap(),
            Ploidy::try_new(2).unwrap(),
        ];
        ploidies.sort();
        assert_eq!(
            ploidies.iter().map(|p| p.get()).collect::<Vec<_>>(),
            vec![1, 2, 4]
        );
    }

    proptest::proptest! {
        /// Over the whole `f64` line: the three rates accept a value **exactly
        /// when** it is a finite number in `[0, 1]`, and an accepted value comes
        /// back bit for bit — no clamp, no round.
        ///
        /// This is what the point assertions above cannot do. Each range check is
        /// two comparisons, and a widened bound on either side is a value reaching
        /// step 8's genotype prior instead of an error — a wrong genotype rather
        /// than a failure. The dense `-2.0..3.0` arm is load-bearing: sampling
        /// `f64::ANY` alone never lands in `(1, 2]`, so it does not notice an
        /// `InbreedingF` whose upper bound has moved to 2.
        #[test]
        fn the_constrained_rates_accept_exactly_the_probabilities_and_round_trip(
            x in proptest::prop_oneof![proptest::num::f64::ANY, -2.0f64..3.0f64]
        ) {
            let is_probability = x.is_finite() && (0.0..=1.0).contains(&x);
            for accepted in [
                ErrorRate::try_new(x).map(ErrorRate::get),
                GenotypeFrequency::try_new(x).map(GenotypeFrequency::get),
                InbreedingF::try_new(x).map(InbreedingF::get),
            ] {
                proptest::prop_assert_eq!(accepted.is_ok(), is_probability, "x = {}", x);
                if let Ok(value) = accepted {
                    proptest::prop_assert_eq!(value.to_bits(), x.to_bits(), "x = {}", x);
                }
            }
        }

        /// Ploidy accepts **every** copy number a genome could have and rejects
        /// only zero. Named `every` and checking three of 255 was the gap: a
        /// ceiling slipped in later would reject a legitimate hexaploid region,
        /// and nothing in the three-value test would notice.
        #[test]
        fn ploidy_accepts_every_non_zero_copy_number_and_round_trips(copies in 0u8..=u8::MAX) {
            match Ploidy::try_new(copies) {
                Ok(ploidy) => {
                    proptest::prop_assert!(copies != 0);
                    proptest::prop_assert_eq!(ploidy.get(), copies);
                }
                Err(rejected) => {
                    proptest::prop_assert_eq!(copies, 0);
                    proptest::prop_assert_eq!(rejected, DomainError::Ploidy(0));
                }
            }
        }

        /// A period is accepted **exactly when** it is in `1..=MAX_MOTIF_LEN`, over a range
        /// reaching well past both the STR scope and a byte — the same lesson `Ploidy`'s
        /// proptest above records. A bound written so that it truncates or masks admits
        /// period 8 while still rejecting 0, 7 and 255, and a test at three points would not
        /// notice; a period no scanner emits would then file loci under a stratum no stutter
        /// model was ever fitted at.
        #[test]
        fn ssr_period_accepts_exactly_the_str_scope(bases in 0usize..=1_000) {
            let inside = (1..=MAX_MOTIF_LEN).contains(&bases);
            proptest::prop_assert_eq!(SsrPeriod::try_new(bases).is_ok(), inside, "bases = {}", bases);
            if let Ok(period) = SsrPeriod::try_new(bases) {
                proptest::prop_assert_eq!(usize::from(period.get()), bases);
            }
        }
    }

    fn genome_position(contig: u32, position: u64) -> GenomePosition {
        GenomePosition {
            contig: ContigId(contig),
            position: Position(position),
        }
    }

    /// The whole reason `GenomePosition` is a struct with this field order and
    /// not a `(Position, ContigId)` pair: sorting must give **genome order** —
    /// contig-major, position-minor — with no comparator written at the call
    /// site. A transposed field order would still compile and still sort; it
    /// would just interleave contigs, which is the failure this test exists to
    /// catch.
    #[test]
    fn sorting_yields_contig_major_position_minor_order() {
        let mut shuffled = vec![
            genome_position(2, 5),
            genome_position(0, 900),
            genome_position(1, 1),
            genome_position(0, 7),
            genome_position(2, 1),
            genome_position(1, 1000),
            genome_position(0, 100),
        ];
        shuffled.sort();

        assert_eq!(
            shuffled,
            vec![
                genome_position(0, 7),
                genome_position(0, 100),
                genome_position(0, 900),
                genome_position(1, 1),
                genome_position(1, 1000),
                genome_position(2, 1),
                genome_position(2, 5),
            ],
            "contig index dominates: contig 2 position 1 sorts after contig 1 position 1000"
        );
    }

    /// Equal keys are legal and compare equal — the order guard rejects only a
    /// strict decrease (`spec/alignment_file.md` §3.2), so `Ord` must report
    /// `Equal` for a repeated position rather than imposing a tie-break of its
    /// own, and `Less` for any genuine advance along the genome.
    #[test]
    fn ordering_is_equal_for_a_repeated_position_and_strict_otherwise() {
        use std::cmp::Ordering;

        let read_start = genome_position(3, 42);
        assert_eq!(read_start.cmp(&genome_position(3, 42)), Ordering::Equal);
        assert!(
            read_start < genome_position(3, 43),
            "later position on same contig"
        );
        assert!(
            read_start < genome_position(4, 1),
            "later contig, earlier position"
        );
    }

    fn region(start: u64, end: u64) -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(0),
            start: Position(start),
            end: Position(end),
        }
    }

    /// **Inclusive is the whole point**, so `len` is where the `+ 1` lives — and
    /// the single-base case is the one that catches a half-open slip.
    #[test]
    fn region_len_counts_both_ends() {
        assert_eq!(region(1, 1).len(), 1, "a single base is length 1, not 0");
        assert_eq!(region(1, 10).len(), 10);
        assert_eq!(
            region(6, 21).len(),
            16,
            "the classification fixture's tract"
        );
        assert!(!region(1, 1).is_empty());
    }

    /// A region can be built inverted (public fields, no constructor), so the
    /// accessors must answer rather than panic.
    #[test]
    fn an_inverted_region_is_empty_not_a_panic() {
        assert!(region(10, 9).is_empty());
        assert_eq!(region(10, 9).len(), 0);
        // Saturating, not wrapping: a wildly inverted span is still 0, not ~u64::MAX.
        assert_eq!(region(u64::MAX, 1).len(), 0);
    }

    #[test]
    fn region_contains_includes_both_bounds() {
        let r = region(6, 21);
        assert!(r.contains(Position(6)), "start is inside");
        assert!(r.contains(Position(21)), "end is inside — inclusive");
        assert!(!r.contains(Position(5)));
        assert!(!r.contains(Position(22)));
    }

    /// Position 0 is representable but meaningless, since 1-based coordinates
    /// begin at one. Recorded rather than enforced: `Position` is unconstrained,
    /// as `MapQual` and `Bp` are, and every place where a 0 is a *bug* rejects it
    /// itself — `RefSeq::fetch_into` with `InvalidStart`, `Locus::new` with
    /// `BadCoordinates`, `classify` with an assert.
    #[test]
    fn position_zero_is_representable_and_rejected_where_it_matters() {
        assert_eq!(Position(0).get(), 0);
        assert!(region(0, 5).contains(Position(0)));
    }

    /// Read-group ids are minted in table order, so their derived `Ord` *is* that
    /// order. Anything reporting a set of read groups can therefore sort them and
    /// get the order the table was built in, which is the order the input files
    /// were given.
    #[test]
    fn read_group_ids_sort_in_table_order() {
        let mut ids = vec![ReadGroupId(2), ReadGroupId(0), ReadGroupId(1)];
        ids.sort();
        assert_eq!(
            ids,
            vec![ReadGroupId(0), ReadGroupId(1), ReadGroupId(2)],
            "the derived Ord is the index order"
        );
    }

    /// Every period the STR scope allows is accepted and reports itself unchanged.
    #[test]
    fn ssr_period_accepts_every_period_in_the_str_scope() {
        for bases in 1..=MAX_MOTIF_LEN {
            let period = SsrPeriod::try_new(bases).expect("inside the STR scope");
            assert_eq!(usize::from(period.get()), bases);
        }
    }

    /// **Zero is the rejection that matters**, because a tract's length becomes a repeat
    /// count by dividing by the period. Seven is the other end of the same range: a
    /// heptamer is outside the scope every other part of ng works in, so admitting one
    /// here would mint a period no scanner ever emits and no stutter model was fitted at.
    /// 258 is the value a `u8` parameter would have silently accepted as a dinucleotide.
    #[test]
    fn ssr_period_rejects_zero_and_anything_past_the_str_scope() {
        assert_eq!(SsrPeriod::try_new(0), Err(DomainError::SsrPeriod(0)));
        assert_eq!(SsrPeriod::try_new(7), Err(DomainError::SsrPeriod(7)));
        assert_eq!(SsrPeriod::try_new(255), Err(DomainError::SsrPeriod(255)));
        assert_eq!(SsrPeriod::try_new(258), Err(DomainError::SsrPeriod(258)));
    }

    /// The two accessors answer the same question in two types, and a motif is the only
    /// thing that mints a period in the live path — so if these ever disagree, every
    /// stratum a locus is filed under moves.
    #[test]
    fn a_motifs_two_period_accessors_agree_at_every_length() {
        for bases in [b"A".as_slice(), b"CA", b"CAG", b"AAAT", b"CACAG", b"ACGTGC"] {
            let motif = Motif::new(bases).expect("a motif inside the STR scope");
            assert_eq!(motif.period(), bases.len());
            assert_eq!(usize::from(motif.ssr_period().get()), motif.period());
        }
    }

    /// A period renders as the bare number, so a message can supply its own word for it.
    #[test]
    fn a_period_renders_as_the_number_of_bases_alone() {
        assert_eq!(
            SsrPeriod::try_new(3).unwrap().to_string(),
            "3",
            "no type name, no unit — the message says 'period {{period}}'"
        );
    }
}
