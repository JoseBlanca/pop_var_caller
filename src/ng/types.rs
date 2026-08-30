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
    ///
    /// **Wrong at the coordinate ceiling, and the saturation above does not cover it:**
    /// `end == u64::MAX` overflows the `+ 1` before anything is subtracted, which is a
    /// panic in a debug build and a length of 0 in the release profile, where overflow
    /// checks are off — a region at the ceiling reporting itself empty. Not reachable
    /// from real contig coordinates, which is why it is recorded rather than fixed here;
    /// it is why [`SampleLocusObservations::reach`] does not obtain its span from this
    /// method, and `locus_generation::tests::a_locus_at_the_coordinate_ceiling_reaches_its_own_end`
    /// is what pins that.
    ///
    /// [`SampleLocusObservations::reach`]: crate::ng::locus_generation::SampleLocusObservations::reach
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

/// Which set of read groups ran together — an index into the run's declared batching.
///
/// A **sequencing batch** is the group of libraries that were sequenced beside one another,
/// as the run was *told*: a flowcell, a plate, a submission. It matters because **a
/// contaminating read is far likelier to come from a neighbour on the same run than from a
/// random member of the species**, so the population a contaminant's genotype is drawn against
/// is the batch and not the cohort (`doc/devel/ng/spec/read_likelihoods.md` §3.6,
/// `doc/devel/ng/arch/parameter_prepass_joint_fit.md` §1.6).
///
/// **It is stated, never inferred.** The grouping is absent from both benchmark cohorts'
/// alignments — the tomato archive's `@RG` lines carry no platform unit, and SRA rewrote the
/// read names — and a pipeline that guessed it from what survives would be wrong in silence.
///
/// **The default is one batch holding the whole run**, which is `BatchId(0)` for every read
/// group. So a run that declares no batching gets the cohort frequency and loses nothing it
/// had, and no consumer branches on the batching's absence.
///
/// Unconstrained for [`ReadGroupId`]'s reason: it indexes a table, an out-of-range value is
/// caught at lookup, and `u32` rather than `u64` because it is an index and not a position.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct BatchId(pub u32);

impl BatchId {
    /// The batch every read group is in when a run declares no batching — **the name the
    /// architecture already uses for it**, `SequencingBatches::all_together`
    /// (`doc/devel/ng/arch/parameter_prepass_joint_fit.md` §1.6).
    pub const ALL_TOGETHER: Self = Self(0);

    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

/// A run's batching **keyed by read group** — entry *i* is [`ReadGroupId`] *i*'s batch.
///
/// **A wrapper rather than a bare `&[BatchId]`, because the sample-keyed batching is the same
/// slice type and means something else.** The two agree in length whenever a run has one
/// library per sample — which is every sample of every benchmark cohort here — so transposing
/// them passes both shape checks and comes back as a wrong contaminant frequency rather than a
/// panic. Sample order and read-group order are minted by different rules, so the mis-key is
/// only invisible, never harmless: the frequency it produces is worth up to 12 nats a read
/// between batches.
///
/// **The same argument the allele-copy views already won here.**
/// `CohortAlleleCopies` and `SampleAlleleCopies` are two types for one shape for exactly this
/// reason, and the measurement recorded there is that the flat-slice version, swapped, silently
/// returned the bare seed at every allele with nothing raised.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct BatchOfEachReadGroup<'a>(pub &'a [BatchId]);

/// A run's batching **keyed by sample** — entry *i* is sample *i*'s batch, in the run's sample
/// order. [`BatchOfEachReadGroup`] says why this is a wrapper and not a slice.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct BatchOfEachSample<'a>(pub &'a [BatchId]);

/// Just the index, for [`ReadGroupId`]'s reason: a message naming a batch supplies its own
/// word for it.
impl fmt::Display for BatchId {
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

/// **How much error the reads behind one observation carry, added up** — the sum over those
/// reads of the natural logarithm of the probability that each is wrong. Always zero or
/// negative, since it sums logarithms of probabilities; production calls it `q_sum`.
///
/// **Held as a count of steps of 1/4,096 of a natural log, not as a float** (the owner,
/// 2026-08-25; `spec/psp_record_encoding.md` §5.1.1). The step is a property of this type,
/// and every route that carries the quantity inherits it — a psp stores the integer it was
/// handed and records the step so a reader can interpret it, and cannot write a file with a
/// step the type did not produce.
///
/// **Why an integer, when the quantity is continuous.** The caller has two routes to the same
/// answer: reading a sample's observations straight from memory, and reading them back from a
/// psp. Those two must produce the same VCF, and that identity is the oracle everything else
/// is checked against. A float stored to full precision breaks it in a way no tolerance
/// repairs, because two routes that add the same terms **in a different order** differ in the
/// last bits — so the check would degrade from *identical* to *within a tolerance*, which is a
/// far weaker test and one that can pass while a chain-id list is being corrupted. Rounding to
/// a step **absorbs** that difference instead: both routes land on the same integer.
///
/// **The precision is not the thing being traded away.** This term goes straight into a
/// genotype likelihood, and the error a step introduces there is the whole risk: a step of
/// 1/16 of a natural log is a 6 % error in that term, 1/256 is 0.4 %, and 1/4,096 is
/// **0.024 %**. The owner took the precision because it costs about 5 % of the file and being
/// wrong about a likelihood costs a genotype.
///
/// **It also deletes a state that used to be reachable.** A not-a-number error sum once came
/// back through `f64::max` as the most confident read the model can express — a confident
/// wrong answer with nothing failing. An integer has no such value.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct SummedLogError(i64);

impl SummedLogError {
    /// Steps to one natural log — the step is `1 / STEPS_PER_NAT`.
    ///
    /// A power of two so that the conversion is exact in both directions for any value a
    /// float can hold exactly at this scale.
    pub const STEPS_PER_NAT: i64 = 4_096;

    /// No error mass at all: the identity for adding these together, and what an observation
    /// with no reads carries.
    pub const NONE: Self = Self(0);

    /// The nearest whole step to `nats`.
    ///
    /// **Round once, at the end of a sum, not once per read.** Rounding each read's own term
    /// and adding the integers would be exactly order-independent, but it accumulates up to
    /// half a step per read — about 0.0012 natural logs at 300 reads a position, ten times the
    /// error of rounding the finished sum. The ordering difference this type exists to absorb
    /// is in the last bits of an `f64`, far below one step, so summing in floating point and
    /// rounding once is both more accurate and enough.
    ///
    /// A value too large for the step count saturates rather than wrapping, and a
    /// not-a-number sum becomes [`NONE`](Self::NONE) — neither is reachable from reads whose
    /// error probabilities are probabilities, and both are stated so that no input produces a
    /// silently wrong number.
    #[inline]
    pub fn from_nats(nats: f64) -> Self {
        if nats.is_nan() {
            return Self::NONE;
        }
        let steps = (nats * Self::STEPS_PER_NAT as f64).round();
        // `as` saturates at the integer's bounds for floats, including infinities.
        Self(steps as i64)
    }

    /// The value in natural logs, for the arithmetic that needs a float.
    #[inline]
    pub fn nats(self) -> f64 {
        self.0 as f64 / Self::STEPS_PER_NAT as f64
    }

    /// The raw step count — what a psp stores and what a header's step interprets.
    #[inline]
    pub fn steps(self) -> i64 {
        self.0
    }

    /// A step count read back from a file, unchanged.
    #[inline]
    pub fn from_steps(steps: i64) -> Self {
        Self(steps)
    }
}

/// Adding two of these is **exact and order-independent**, which is the property the whole
/// type exists for: a cohort merge folding one read's evidence from several records, or a
/// selection pooling the alleles it cut into a leftover, gets the same answer whatever order
/// it works in. Saturating rather than wrapping, for the reason [`Self::from_nats`] gives.
impl std::ops::Add for SummedLogError {
    type Output = Self;

    #[inline]
    fn add(self, other: Self) -> Self {
        Self(self.0.saturating_add(other.0))
    }
}

impl std::ops::AddAssign for SummedLogError {
    #[inline]
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl std::iter::Sum for SummedLogError {
    #[inline]
    fn sum<I: Iterator<Item = Self>>(terms: I) -> Self {
        terms.fold(Self::NONE, |running, term| running + term)
    }
}

/// `-3.5 ln` — the unit is in the rendering, because a bare number here reads as a
/// probability, a Phred score or a step count depending on who is looking.
impl fmt::Display for SummedLogError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} ln", self.nats())
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

/// Which allele of one locus: an index into that locus's candidate-allele
/// table (`CandidateAlleles`, `doc/devel/ng/arch/calling_em_loop.md` §2), where
/// index `0` is the reference allele and is always present.
///
/// **An id means nothing without the locus it was minted at.** Allele 1 here
/// and allele 1 at the next locus are different pieces of sequence, exactly as
/// a [`Position`] names no base until a [`ContigId`] is known. The type carries
/// no locus, so an id must not outlive the table it indexes — the calling loop
/// keeps the two together and mints the owned `Genotype` multiset, one id per
/// genome copy, only at the end (`arch/calling_em_loop.md` §2).
///
/// Unconstrained — any `u16` is a legal index at the type level, and an
/// out-of-range id is caught when the table is read — so the field is public
/// and there is no checked constructor. `u16` and not `u32`, because a locus is
/// pruned to a handful of candidates: production keeps 6 alleles per record by
/// default and refuses to be configured above 16
/// (`DEFAULT_MAX_ALLELES_PER_RECORD` and `MAX_ALLELES_PER_VAR_CAP`,
/// `var_calling::per_group_merger`), so the ceiling here is about four thousand
/// times the widest cap that can be asked for.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct AlleleId(pub u16);

impl AlleleId {
    /// The reference allele — index `0` of every locus's candidate table, always
    /// present (`arch/calling_em_loop.md` §2). Named so the convention is
    /// greppable and no consumer spells it as a bare `0`: the reference is what
    /// every downstream branch tests against, REF against ALT in the VCF and the
    /// homozygous-reference genotype in the prior.
    pub const REFERENCE: Self = Self(0);

    /// Whether this id names the reference allele.
    #[inline]
    pub fn is_reference(self) -> bool {
        self == Self::REFERENCE
    }

    #[inline]
    pub fn get(self) -> u16 {
        self.0
    }
}

/// A quality on the Phred scale — `-10 log10(p)`, where `p` is the chance that
/// the thing being scored is wrong. 20 means one call in a hundred is wrong,
/// 30 one in a thousand, 60 one in a million.
///
/// **This is the scale ng writes, not the one it works in.** The internal
/// currency is the natural logarithm ([`LogProb`]); Phred exists because VCF's
/// `QUAL` and `GQ` columns are written on it. The two are logarithms to
/// different bases with opposite signs, so one type for both is how a
/// log-probability ends up added to a quality and read as a plausible wrong
/// number instead of failing to compile. Every crossing between the two scales
/// is a named function ([`Self::from_log_prob`]) and never a bare `as` cast,
/// which would change the width of a number while saying nothing about its
/// scale. (The narrowing to `f32` inside that function is a width change and
/// nothing else — the scaling is the multiply that precedes it.)
///
/// Constrained, so validated: the field is private and construction goes
/// through [`Self::try_new`]. Below zero is a probability above one, and a
/// caller's arithmetic gone wrong. Infinite is a probability of exactly zero —
/// which log space represents happily (`ln 0 = -∞`, the score of an impossible
/// read line-up: see [`LogProb`]) and this scale has no number for at all. The
/// two are **not** the same kind of event, so they do not share an error: the
/// first is [`DomainError::Phred`], the second [`DomainError::PhredInfinite`].
///
/// **Neither is clamped here, because where to cap is the consumer's call.**
/// Production caps `GQ` at
/// [`DEFAULT_MAX_GQ_PHRED`](crate::var_calling::posterior_engine::DEFAULT_MAX_GQ_PHRED)
/// — 99, the GATK and bcftools convention, configurable up to
/// [`GQ_PHRED_RANGE_MAX`](crate::var_calling::posterior_engine::GQ_PHRED_RANGE_MAX)
/// — at the point it fills the column, and pins the posterior just below one
/// first so the infinity rarely arises. A clamp inside the type would pick a
/// ceiling for every future consumer and hide the arithmetic that produced the
/// value.
///
/// **A quality ng computed itself has no constructor here yet.**
/// `arch/ng_step_interfaces.md` §1 allows one for that source — a `new` that
/// `debug_assert!`s the bound and clamps only a float-epsilon overrun — and the
/// step that first fills a `GQ` column is where it should land, rather than as a
/// clamp written at that call site.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct Phred(f32);

impl Phred {
    /// The one check, and every constructor goes through it. A quality below
    /// zero or a `NaN` says the caller's arithmetic went wrong; an infinite one
    /// says the scored probability is exactly zero, which is a different event
    /// and gets [`DomainError::PhredInfinite`]. Nothing is coerced either way.
    ///
    /// `NaN` needs no check of its own: no comparison with it is true, so
    /// `quality >= 0.0` already rejects it, the way `MismatchFraction`'s range
    /// check does.
    ///
    /// Zero has one spelling. `-0.0` passes `>= 0.0` — it *is* zero — but it
    /// prints as `-0`, and [`Self::from_log_prob`] produces exactly it at
    /// `ln 1`, since `-k * 0.0` is `-0.0`. A `QUAL` column cannot carry a minus
    /// sign on a certainty, so the sign is normalised here, at the one door in.
    pub fn try_new(quality: f32) -> Result<Self, DomainError> {
        if quality == f32::INFINITY {
            return Err(DomainError::PhredInfinite);
        }
        (quality >= 0.0 && quality.is_finite())
            .then_some(Self(if quality == 0.0 { 0.0 } else { quality }))
            .ok_or(DomainError::Phred(quality))
    }

    /// The named crossing from ng's working scale to the output one: `ln(p)` in,
    /// `-10 log10(p)` out.
    ///
    /// `Err` wherever the Phred scale cannot follow the logarithm, and the two
    /// causes are told apart by the error rather than by the caller inspecting a
    /// float:
    ///
    /// - `ln(p) = -∞` — probability zero, an infinite quality.
    ///   [`DomainError::PhredInfinite`]. Not a bug: `-∞` is inside [`LogProb`]'s
    ///   documented domain, and a consumer's answer is to cap at its own ceiling
    ///   and carry on.
    /// - `ln(p) > 0` — a probability above one, so a negative quality.
    ///   [`DomainError::Phred`], carrying the quality that was computed rather
    ///   than the log probability handed in; divide by `10/ln(10)` to recover it.
    ///
    /// The scaling is done at `f64`, the width [`LogProb`] holds, and narrowed
    /// once at the end; narrowing first would round the log probability before
    /// scaling it. The narrowing **saturates** rather than wrapping, so a log
    /// probability finite in `f64` but below about `-7.8e37` becomes `+∞` and is
    /// refused as an infinite quality — the same answer `ln(p) = -∞` gets, for
    /// the same reason: the scale has no number for it. Unreachable from real
    /// data, and stated so the `is_finite` guard is not simplified away.
    pub fn from_log_prob(log_p: LogProb) -> Result<Self, DomainError> {
        Self::try_new((-PHRED_PER_NAT * log_p.get()) as f32)
    }

    #[inline]
    pub fn get(self) -> f32 {
        self.0
    }
}

/// Phred per nat — per unit of natural logarithm: `10 / ln(10)`, which is
/// `10 log10(e)`. The whole of the [`LogProb`] → [`Phred`] conversion, written
/// once. A log probability is negative, so the crossing negates as it scales
/// (`-PHRED_PER_NAT * ln(p)`), which is how htslib and this repository's BAQ
/// port both spell it.
///
/// **A second constant of this name exists** — `baq::probaln::PHRED_PER_NAT`,
/// the same quantity as the four-digit literal htslib compiled, `4.343`. That
/// one is kept at htslib's value on purpose, because the BAQ port has to
/// reproduce htslib's numbers byte for byte; this one is the full ratio,
/// because nothing here is reproducing another program's arithmetic. Do not
/// unify them.
const PHRED_PER_NAT: f64 = 10.0 / std::f64::consts::LN_10;

// ---------------------------------------------------------------------
// The parameters a caller runs on — four scalars step 4 fits (the error
// rate, the genotype frequencies, the inbreeding coefficient, the expected
// heterozygosity) and one it is handed (ploidy). Five types and not one
// shared `Probability`: four of them are fractions — three closed at
// `[0, 1]`, the inbreeding coefficient half-open at `[0, 1)` — and no range
// tells them apart in a way a compiler could use, so a single type would let
// an inbreeding coefficient be handed to something expecting an error rate
// and compile (`arch/parameter_prepass_generic.md` §2.1).
//
// They live in the shared vocabulary rather than in `parameter_estimation/`
// because their consumers are *other* steps: the likelihood (step 7) reads
// the error rate, the genotype prior (step 8) reads the genotype
// frequencies, the inbreeding coefficient and the expected heterozygosity,
// and ploidy reaches both. **Step 7 has no module yet and step 8's holds no
// code** — defining the types here is what keeps them from importing out of a
// sibling step's module when they arrive. They are also the seed of the
// `genotype`/`params` split `module_layout.md` principle 3 anticipates for
// this file.
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
/// (`parameter_estimation::ssr::slippage`), so **seven** constructors now share this
/// predicate. It is not hypothetical drift — `SiteNoise::try_new`
/// (`parameter_estimation::generic`) already spells the range test by hand, which is the
/// eighth probability in the crate and the one this predicate did not reach.
///
/// **[`InbreedingF`] is `[0, 1)` and is not an instance of that drift**: it composes this
/// predicate and rejects the ceiling on top of it, so the closed range is still written
/// once and its one exception is written where the exception is explained.
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
///
/// **Half-open: `[0, 1)`. What that buys is finiteness, not a small number.**
/// The genotype prior multiplies its heterozygote branch by `1 − F`, so at
/// `F = 1` the factor is exactly zero and `ln(1 − F)` is `−∞`: every heterozygote
/// becomes impossible and no read evidence, however clean, could produce a
/// heterozygous call. Excluding the endpoint is what keeps that branch a finite
/// number. It is **not** a numerically meaningful cap — the largest value the
/// type accepts still leaves `1 − F = 2⁻⁵³`, which is about 160 on the Phred
/// scale against every heterozygote, where two clean alternative bases at Q30
/// supply 60 (`spec/calling_priors.md` §7).
///
/// **Capping the estimate is a different job and belongs to whoever fits one.**
/// Production's estimator clamps its *fitted* value at `0.99`
/// ([`MAX_INBREEDING_COEFFICIENT`](crate::paralog::inbreeding::MAX_INBREEDING_COEFFICIENT)) —
/// 20 Phred, which read evidence can overcome. Its `--inbreeding-coefficient`
/// flag is a second door and is not clamped: the parser accepts the closed
/// `[0, 1]` (`pop_var_caller::cli::parsers::parse_inbreeding_coefficient`) and
/// the value goes to the engine as given (`var_calling::pipeline`). That is the
/// gap this newtype closes for ng — the limit is unrepresentable here whichever
/// door a value came through.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct InbreedingF(f64);

impl InbreedingF {
    /// The only constructor. A coefficient that is not a fraction in `[0, 1)` is
    /// rejected rather than coerced.
    ///
    /// **Two rejections, not one**, because they mean different things to whoever
    /// reads the message: a value outside `[0, 1]` is not a fraction at all
    /// ([`DomainError::InbreedingF`]), while exactly `1` is a fraction that this
    /// type deliberately excludes ([`DomainError::InbreedingFAtCeiling`]). The
    /// shared `[0, 1]` predicate is reused for the first and the ceiling is checked
    /// on top of it: [`checked_probability`] admits `1.0` for [`ErrorRate`] and
    /// [`GenotypeFrequency`], where it is a real answer, and tightening it there
    /// would break both (`arch/calling_priors.md` §2.1).
    pub fn try_new(coefficient: f64) -> Result<Self, DomainError> {
        let fraction = checked_probability(coefficient, DomainError::InbreedingF)?;
        if fraction == 1.0 {
            return Err(DomainError::InbreedingFAtCeiling(fraction));
        }
        Ok(Self(fraction))
    }

    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// How different two chromosomes drawn at random from the cohort are expected to be
/// at an ordinary site — the chance they differ there, averaged over the sites this
/// caller treats as ordinary. The genotype prior's `θ` (`spec/calling_priors.md` §4).
/// A probability in `[0, 1]`.
///
/// **Differ, not "carry different bases"**: one `θ` covers substitutions and short
/// insertions and deletions alike, because the pre-pass measures one number for both
/// (`spec/calling_priors.md` §4.2).
///
/// **Not the non-reference rate**, which counts how often the cohort differs from the
/// *reference sequence* and so books every quirk of the one accession the reference
/// was assembled from as cohort polymorphism. The two are different numbers on any
/// panel whose reference is not one of its own members, and a prior built from the
/// second would claim a diversity the population does not have
/// (`spec/calling_priors.md` §4).
///
/// **Ordinary sites only.** Repeat tracts mutate orders of magnitude faster, so how
/// variable they are is a separate measurement and this value must never stand in for
/// it (`spec/calling_priors.md` §5). Production's STR path never measured it at all: it
/// hardcodes freebayes' default `SFS_THETA = 0.01`
/// (`src/ssr/cohort/freebayes_emit.rs`), marked there "Fixed, not a per-run knob" — and
/// that number is a population-scaled mutation rate, a different quantity in different
/// units from a heterozygosity in `[0, 1]`. Not repeating that is why the two are
/// separate types rather than one shared float.
///
/// Source: the pre-pass, which has two routes and today supplies this from one of them.
/// The joint fit reads it off its fitted density (`JointFit::expected_heterozygosity`,
/// `parameter_estimation::joint`) and runs at every cohort size down to one sample. The
/// per-sample histogram route (`parameter_estimation::generic`) supplies the ingredient
/// rather than the number — each sample's *observed* heterozygosity, of which `θ` is the
/// mean of `Hobs / (1 − F)` across samples (`spec/calling_priors.md` §4) — and **nothing
/// computes that mean yet**. Where no fit exists at all — too few sites, or no `F` for
/// the sample — the caller falls back to [`Self::SPECIES_FALLBACK`] and must say so in
/// its output.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct ExpectedHeterozygosity(f64);

impl ExpectedHeterozygosity {
    /// What a run assumes when nothing could be fitted: roughly human nucleotide
    /// diversity, one difference per thousand bases.
    ///
    /// **The value of last resort, and a run that lands on it must say so in its
    /// output** — a run on a species-range guess and a run on a measured diversity are
    /// otherwise indistinguishable (`spec/calling_priors.md` §4). The thing that will
    /// carry that is the genotype prior's `SeedRegime::FallbackDiversity`
    /// (`arch/calling_priors.md` §2.3), and **any code path that reads this constant
    /// owes it**. Nothing reads it yet.
    ///
    /// **It must be overridable, and no door exists yet.** No command-line flag or
    /// configuration field sets it; the run configuration that will is the calling
    /// loop's (`doc/devel/ng/impl_plan/calling_loop.md`). Until then a run on an
    /// unusual species has no way to correct it.
    ///
    /// **It is a human figure, and it is not a starting point for another species.**
    /// Which way it is wrong is not fixed: this project's own tomato panel fits *below*
    /// it — 6 differences per 10,000 bases, against the 10 per 10,000 here
    /// (`spec/calling_priors.md` §4.1) — while a diverse outcrosser would sit above.
    ///
    /// Taken from production's `DEFAULT_DIVERSITY_PRIOR`
    /// (`src/var_calling/diversity.rs`), value and reasoning. Production is frozen and
    /// this constant is ng's own: the two are not tied, and ng may move this one.
    ///
    /// **A value of the type rather than a bare `f64`**, following [`AlleleId::REFERENCE`]
    /// — the same reason. As a loose float it constructs an [`ErrorRate`], an
    /// [`InbreedingF`] and a [`GenotypeFrequency`] just as happily, which is the
    /// confusion the five separate types in this section exist to prevent, and it would
    /// be the only diversity constant in the shared vocabulary for an STR-path author to
    /// reach for.
    pub const SPECIES_FALLBACK: Self = Self(1e-3);

    /// The only constructor. A heterozygosity that is not a probability in `[0, 1]` is
    /// rejected rather than coerced.
    pub fn try_new(heterozygosity: f64) -> Result<Self, DomainError> {
        checked_probability(heterozygosity, DomainError::ExpectedHeterozygosity).map(Self)
    }

    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// **How often a chromosome drawn at random from the population carries something other
/// than the reference base**, at an ordinary site — the mean alternative-allele
/// frequency. A probability in `[0, 1]`.
///
/// **The partner of [`ExpectedHeterozygosity`], and a type of its own so the two cannot
/// be swapped.** The SNP/indel genotype prior's two concentration numbers are exactly
/// this frequency and a total conviction in other clothes, and
/// `doc/devel/ng/spec/ordinary_site_seed.md` §3's identity turns the pair back into
/// them. Both are probabilities in `[0, 1]`, both are population quantities with no
/// panel in them, and both reach the seed builder in the same call — so as bare floats
/// a swapped pair compiles and returns a seed that is wrong in a way no downstream check
/// refuses.
///
/// **They are different questions about the same population.** This one asks how often
/// *one* chromosome is non-reference; the heterozygosity asks how often *two* differ.
/// On a population where the alternative allele is rare the second is about twice the
/// first, and on one where the reference base is the rare one it is far smaller than
/// the first — so their sizes do not separate them either.
///
/// Source: the joint fit's own fitted density, in closed form
/// (`FrequencyDensity::expected_alternative_frequency`,
/// `doc/devel/ng/spec/ordinary_site_prior_moments.md` §2). A run whose pre-pass fitted no
/// density has none of this, and its seed falls to the neutral shape at whatever
/// diversity it did fit.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct ExpectedAlternativeFrequency(f64);

impl ExpectedAlternativeFrequency {
    /// The only constructor. A frequency that is not a probability in `[0, 1]` is
    /// rejected rather than coerced.
    pub fn try_new(frequency: f64) -> Result<Self, DomainError> {
        checked_probability(frequency, DomainError::ExpectedAlternativeFrequency).map(Self)
    }

    #[inline]
    pub fn get(self) -> f64 {
        self.0
    }
}

/// The chance that two repeat-tract copies drawn at random from the cohort carry
/// different numbers of repeats — Nei's gene diversity, measured on repeat tracts.
///
/// **The same question as [`ExpectedHeterozygosity`], asked of a different site class,
/// and it is a type of its own so the two cannot be swapped.** Repeat tracts mutate
/// orders of magnitude faster than bases do, so the two numbers are not the same size
/// and neither substitutes for the other; the pre-pass is required to emit them
/// separately precisely so a consumer cannot confuse them
/// (`doc/devel/ng/spec/calling_priors.md` §5,
/// `doc/devel/ng/spec/parameter_prepass_cohort.md` §3). Today's production STR path is
/// what happens when they are not separated: it takes a fixed `SFS_THETA = 0.01`
/// (`src/ssr/cohort/freebayes_emit.rs`), freebayes' SNP-scale default, commented
/// *"Fixed, not a per-run knob"* — a constant standing in for a quantity nobody
/// measured. As bare floats, that substitution compiles.
///
/// **A gene diversity and not a heterozygosity, because a tract is multi-allelic by
/// length.** "How often do two copies differ" is still the right question; the
/// two-allele arithmetic is not. Like [`ExpectedHeterozygosity`] it is a cohort
/// quantity, so an individual's inbreeding does not enter it — the pre-pass divides the
/// observed rate by `(1 − F)` before it gets here.
///
/// **Fitted at every cohort size down to one.** With a single genome it is that
/// genome's own observed rate over its `(1 − F)`; one diploid genome carries two copies
/// of every tract, and how often those two differ is exactly what this asks. So unlike
/// a frequency spectrum it never comes back absent for want of a panel.
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug)]
pub struct RepeatGeneDiversity(f64);

impl RepeatGeneDiversity {
    /// The only constructor. A value that is not a probability in `[0, 1]` is rejected
    /// rather than coerced.
    ///
    /// **The whole of `[0, 1]` is accepted**, and there is no consumer left that can
    /// refuse a value inside it. Until 2026-08-26 the repeat tract's genotype prior
    /// scaled a constructed geometric shape to reproduce this measurement, which was
    /// possible only below a ceiling the shape itself set — at one outbred sample, below
    /// every measurement — and the builder refused above it. The prior is now seeded from
    /// the length spectrum the joint repeat fit produces per stratum
    /// (`doc/devel/ng/spec/population_diversity.md` §4.2), which asserts no such scaling,
    /// so nothing reads this value today.
    ///
    /// **It stays because the pre-pass still owes it**
    /// (`doc/devel/ng/spec/parameter_prepass_cohort.md` §3): the STR gene diversity is one
    /// of the two diversities that step is specified to emit, separately from
    /// [`ExpectedHeterozygosity`] and precisely so that a consumer cannot confuse them.
    pub fn try_new(diversity: f64) -> Result<Self, DomainError> {
        checked_probability(diversity, DomainError::RepeatGeneDiversity).map(Self)
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
/// first variant; later constrained types add their own variants as they arrive.
/// `#[non_exhaustive]` so matchers accept those future variants without breaking.
///
/// **`PartialEq` is IEEE equality on the float payloads, so an error carrying a
/// `NaN` is not equal to itself.** A `NaN` input is exactly what the four rate
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
    /// An [`InbreedingF`] was built from a value outside `[0, 1)` — and not from
    /// the ceiling itself, which gets [`Self::InbreedingFAtCeiling`].
    ///
    /// **The message names `[0, 1)`, the type's whole range, not the part this
    /// variant covers.** Someone who mistypes `1.5` reads the range and retries;
    /// if it said `[0, 1]` they would retry with `1.0` and be refused again. Every
    /// value this variant can carry is outside `[0, 1)` too, so nothing is lost by
    /// naming the wider truth.
    #[error("inbreeding coefficient {0} is not a finite fraction in [0, 1)")]
    InbreedingF(f64),
    /// An [`InbreedingF`] was built from exactly `1` — a fraction, but the one
    /// value the half-open type excludes.
    ///
    /// **Its own variant rather than one message for both**, because the two
    /// rejections need different remedies. `1.5` is a typo and naming the range is
    /// the whole answer. `1.0` is a coherent request the model refuses, and the
    /// answer is what refusing it means: a prior under which no sample can ever be
    /// called heterozygous. Carrying that explanation on the shared variant would
    /// put it in front of someone who typed `-0.5`, where it is noise.
    #[error(
        "inbreeding coefficient {0} would make every heterozygote impossible; \
         the accepted range is [0, 1), and an estimate should sit well below it \
         (production caps a fitted one at 0.99)"
    )]
    InbreedingFAtCeiling(f64),
    /// An [`ExpectedHeterozygosity`] was built from a value that is not a finite
    /// probability in `[0, 1]`.
    ///
    /// Its own variant beside [`Self::GenotypeFrequency`], and here the two are
    /// genuinely easy to confuse, because **ng carries a heterozygosity under each
    /// type**. This one draws its two chromosomes from the **cohort**, so an
    /// individual's inbreeding does not touch it; the one a [`GenotypeFrequency`]
    /// carries (`parameter_estimation::generic`'s `observed_heterozygosity`) draws them
    /// from **one individual**, so inbreeding drives it down — the two differ by a
    /// factor of `(1 − F)`. A message naming the wrong one sends the reader to the
    /// wrong fit.
    #[error("expected heterozygosity {0} is not a finite probability in [0, 1]")]
    ExpectedHeterozygosity(f64),
    /// An [`ExpectedAlternativeFrequency`] was built from a value that is not a finite
    /// probability in `[0, 1]`.
    ///
    /// **Its own variant beside [`Self::ExpectedHeterozygosity`], because the two arrive
    /// together.** They are the two numbers the SNP/indel genotype prior is built from
    /// and they reach the seed builder in one call, so a message naming the wrong one
    /// sends the reader to the wrong half of a two-number fit.
    #[error("expected alternative-allele frequency {0} is not a finite probability in [0, 1]")]
    ExpectedAlternativeFrequency(f64),
    /// A [`RepeatGeneDiversity`] was built from a value that is not a finite
    /// probability in `[0, 1]`.
    ///
    /// **Its own variant beside [`Self::ExpectedHeterozygosity`], and the pair is the
    /// point.** Both are "how often do two copies drawn from the cohort differ", and
    /// they are measured on different site classes and come out orders of magnitude
    /// apart — bases mutate slowly, repeat tracts fast. A message naming the wrong one
    /// sends the reader to the wrong half of the pre-pass.
    #[error("repeat gene diversity {0} is not a finite probability in [0, 1]")]
    RepeatGeneDiversity(f64),
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

    /// A locus's reads were laid out across the offset buckets in a number a
    /// locus shape cannot record: none at all, or more than the read cap
    /// (`parameter_estimation::ssr::stratum_table::LocusShape`).
    ///
    /// **Carries `u64` where the shape stores `u8`**, for the reason
    /// [`Self::SsrPeriod`] carries a `usize`: the value the message must name is
    /// the one the *caller* offered, and a caller counting a deep locus's reads
    /// holds a wide integer. A variant as narrow as the storage is what makes
    /// `try_new(counts as u8)` the natural call, under which 260 reads arrive as
    /// 4 and validate.
    ///
    /// **The cap travels in the error rather than being named in the message
    /// text.** Two reasons, and the first alone would not settle it: this file is
    /// the shared vocabulary and the cap is the STR path's own constant, so a
    /// `{MAX_LOCUS_READS}` here would point the wrong way up the module tree.
    /// The second is that the cap is a value the design expects to move —
    /// `arch/parameter_prepass_ssr.md` §7 leaves it as an impl-time confirmation
    /// against the precision it costs — so an error that carries the number in
    /// force says what a run rejected against rather than what this file was
    /// compiled believing.
    #[error("a locus shape holds {reads} reads, outside the 1..={cap} one can record")]
    SsrLocusShapeReads {
        /// How many reads the caller offered, over all buckets and the guard.
        reads: u64,
        /// The cap in force — `MAX_LOCUS_READS`.
        cap: u32,
    },

    /// A locus offered more mismatched bases than bases compared
    /// (`parameter_estimation::ssr::stratum_table::BaseComparison`).
    ///
    /// **Not a noisy locus — a swapped pair.** The two counts have the same type
    /// and one is a subset of the other, so the failure this catches is a caller
    /// handing them over the wrong way round. Left standing, what it does depends
    /// on what the locus is pooled with, and the milder-looking outcome is the
    /// dangerous one: alone it drives the stratum's rate above one, which
    /// [`ErrorRate`] refuses, so the run dies where the rate is read; pooled with
    /// well-formed loci it survives as a plausible wrong rate that nothing
    /// downstream can question.
    #[error(
        "a locus counted {bases_mismatched} mismatched bases among only {bases_compared} \
         compared — mismatched bases cannot outnumber compared ones, so the two counts are \
         swapped or one was miscounted"
    )]
    SsrBaseComparison {
        /// Bases of the locus's reads that were compared against the tract.
        bases_compared: u32,
        /// Of those, the ones the caller says differed.
        bases_mismatched: u32,
    },

    /// A [`Phred`] was built from a negative value or a `NaN` — the caller's
    /// arithmetic went wrong.
    ///
    /// **Its own variant, and the message says "quality" rather than
    /// "probability"**, for the reason the rate variants above have one each: a
    /// Phred is a probability re-expressed as a logarithm, so a message naming a
    /// probability sends the reader hunting for a number between zero and one
    /// when the number in hand is a 30 or a −4.
    #[error("phred quality {0} is not a number at or above zero")]
    Phred(f32),

    /// A [`Phred`] was built from positive infinity — a scored probability of
    /// exactly zero, which the Phred scale has no number for.
    ///
    /// **Its own variant because it is not a caller bug**, and that is the
    /// distinction [`Self::Phred`] cannot carry. `-∞` is inside [`LogProb`]'s
    /// documented domain — the score of an impossible read line-up, which log
    /// space carries on purpose — so a consumer meeting this has a routine
    /// answer: cap at its own ceiling and carry on, as production does at
    /// `DEFAULT_MAX_GQ_PHRED`. Sharing one variant would make every such
    /// consumer tell a routine cap from a broken sum by testing the payload's
    /// sign and finiteness.
    ///
    /// No payload: the value is always `+∞`, so there is nothing to report that
    /// the variant's own name does not already say.
    #[error("phred quality is infinite — the scored probability is exactly zero")]
    PhredInfinite,
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

// ---------------------------------------------------------------------
// The genotype — calling's output vocabulary, shared across steps
// ---------------------------------------------------------------------

/// One individual's genotype at one locus: which alleles that individual
/// carries, one entry per copy of the genome. A **multiset** — order does not
/// matter and repeats are the point, since carrying two copies of the same
/// allele is what homozygous means.
///
/// **This is an output, not a working value.** The calling loop's currency is a
/// row index into the locus's genotype table, because every sample at a locus is
/// scored over the same candidate genotypes and an index is what a score array
/// is addressed by. A `Genotype` is minted from a row only on the loop's last
/// pass, when the locus's calls are written out, which is why this type is small
/// and owns no arithmetic (`arch/calling_em_loop.md` §2).
///
/// The field is **private, and holding the alleles sorted is the reason.** Two
/// genotypes that name the same alleles must compare equal whichever order they
/// were built in, and the cheapest way to get that from derived `PartialEq`,
/// `Ord` and `Hash` is to have one spelling — so [`Self::new`] sorts, and
/// nothing else can construct one. Privacy also keeps ploidy out of the
/// surface: diploid is simply two entries, and a polyploid region changes what
/// the caller passes to [`Self::new`], not this type or anything that consumes
/// it (`arch/ng_step_interfaces.md` §2).
///
/// **The derived [`Ord`] sorts by allele, lowest first, and then by length** —
/// it is `Box<[AlleleId]>`'s lexicographic order over the sorted entries, so
/// `0/0` precedes `0/1` precedes `1/1`, and at mixed ploidy a shorter genotype
/// precedes a longer one sharing its prefix (`[0]` before `[0, 0]`). It exists
/// to give a deterministic output order, not to rank genotypes by anything
/// genetic.
///
/// **How many alleles it holds is the ploidy at that region**, read as
/// `genotype.alleles().len()`. There is no `ploidy()` accessor: one returning a
/// bare `u8` would hand back a number that no longer says what it counts, which
/// is the whole reason [`Ploidy`] is a type rather than an integer — a bare `u8`
/// ploidy is interchangeable at the type level with a bare `u8` mapping quality
/// or base quality, and this file has three such types. One returning a
/// [`Ploidy`] would be no better: `Ploidy` refuses zero copies, so the accessor
/// would have to be fallible for a case [`Self::new`] already refuses outright.
///
/// **And no `is_homozygous()`, deliberately**, though the interfaces sketch had
/// one. `GenotypeTable::homozygous_allele_for` is the *one* homozygous test
/// (`arch/calling_priors.md` §3.2): nothing else may decide homozygosity, so
/// that the rule for above diploidy has a single place to change. A second test
/// here would be the place it silently diverges.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Genotype(Box<[AlleleId]>);

impl Genotype {
    /// The only constructor, and it **sorts** — see the type's doc comment for
    /// why one spelling per genotype is what makes the derived `PartialEq`,
    /// `Ord` and `Hash` mean what a reader expects. The argument is taken by
    /// value, so no caller can observe its own vector reordered.
    ///
    /// # Panics
    ///
    /// On an empty multiset, which is not a genotype at any ploidy: [`Ploidy`]
    /// refuses zero copies because a genome with none is not a genome, and this
    /// type would otherwise make the same quantity legal by the back door. The
    /// check is one comparison per sample per locus, on the last pass only, so
    /// it is not on the loop's hot path — and what it stops is a `GT` field
    /// naming no allele at all, written to a VCF with nothing between here and
    /// the writer having objected.
    ///
    /// Takes an owned `Vec` so the sort happens in place. `into_boxed_slice`
    /// then reuses the buffer when the vector is exactly sized and reallocates
    /// when it is not — a caller that pushed one id per genome copy into a fresh
    /// `Vec` pays one copy of a handful of `u16`s, which is why the signature is
    /// chosen for clarity rather than for that. `sort_unstable` because these
    /// are plain indices: two entries that compare equal are the same bit
    /// pattern, so no fixture could tell a stable sort from an unstable one, and
    /// the stable one would allocate.
    pub fn new(mut alleles: Vec<AlleleId>) -> Self {
        assert!(
            !alleles.is_empty(),
            "a genotype holds one allele per genome copy, and the smallest genome has one \
             copy — an empty multiset is not a haploid call, it is a sample with no genome"
        );
        alleles.sort_unstable();
        Self(alleles.into_boxed_slice())
    }

    /// The alleles carried, in sorted order, one entry per copy of the genome.
    #[inline]
    pub fn alleles(&self) -> &[AlleleId] {
        &self.0
    }
}

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
        assert_eq!(AlleleId(0).get(), 0);
        assert_eq!(AlleleId(u16::MAX).get(), u16::MAX);
    }

    /// Index `0` is the reference allele at every locus, and the constant is
    /// what stops each consumer spelling that convention as a bare `0`.
    #[test]
    fn allele_id_zero_is_the_reference_allele() {
        assert_eq!(AlleleId::REFERENCE.get(), 0);
        assert!(AlleleId::REFERENCE.is_reference());
        assert!(!AlleleId(1).is_reference());
        assert!(!AlleleId(u16::MAX).is_reference());
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

    /// The three closed `[0, 1]` rates accept both endpoints — a genotype frequency
    /// of exactly zero and a fully invariant cohort's expected heterozygosity of
    /// exactly zero are real answers, so a half-open check would reject valid data
    /// there.
    ///
    /// [`InbreedingF`] is the exception: it is `[0, 1)`, so it accepts its lower
    /// endpoint only, and its rejection of the upper one is asserted in
    /// [`each_constrained_rate_rejects_out_of_range_in_both_directions`] below. What
    /// is asserted here is that **exactly one value** is excluded — the very next
    /// `f64` down still constructs, so the type removes the mathematical limit and
    /// nothing else. Keeping an estimate away from that limit is a separate job and
    /// belongs to whoever fits one.
    #[test]
    fn each_constrained_rate_accepts_the_endpoints_of_its_own_range() {
        assert_eq!(ErrorRate::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(ErrorRate::try_new(1.0).unwrap().get(), 1.0);
        assert_eq!(ErrorRate::try_new(0.001).unwrap().get(), 0.001);

        assert_eq!(GenotypeFrequency::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(GenotypeFrequency::try_new(1.0).unwrap().get(), 1.0);

        assert_eq!(ExpectedHeterozygosity::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(ExpectedHeterozygosity::try_new(1.0).unwrap().get(), 1.0);

        assert_eq!(
            ExpectedAlternativeFrequency::try_new(0.0).unwrap().get(),
            0.0
        );
        assert_eq!(
            ExpectedAlternativeFrequency::try_new(1.0).unwrap().get(),
            1.0
        );

        assert_eq!(InbreedingF::try_new(0.0).unwrap().get(), 0.0);
        let just_below_one = f64::from_bits(1.0f64.to_bits() - 1);
        // Pinned, not assumed: any value below one would satisfy the assertion that
        // follows, and only the nearest one shows that a single value is excluded.
        assert_eq!(f64::from_bits(just_below_one.to_bits() + 1), 1.0);
        assert_eq!(
            InbreedingF::try_new(just_below_one).unwrap().get(),
            just_below_one
        );
    }

    /// **Neither of the coefficient's rejections may tell a user the accepted range
    /// is `[0, 1]`.** Someone who mistypes `1.5`, reads that, and retries `1.0` is
    /// refused again — the range in the message is the only thing they have to
    /// retry from. That is the property; how the message words the range is not
    /// pinned, so a rewrite that still refuses to name the closed range is free.
    ///
    /// The separate ceiling variant earns its place on what it *adds*: it says what
    /// refusing `1` means, which is a prior under which no sample can ever be
    /// called heterozygous. The other message must not carry that, because it is
    /// noise in front of someone who typed `-0.5`.
    ///
    /// **Both messages are fetched through [`InbreedingF::try_new`]**, not built by
    /// hand, so this also pins which variant the constructor picks — a constructor
    /// returning the wrong one would otherwise be caught by a single assertion in
    /// one other test.
    #[test]
    fn neither_inbreeding_rejection_names_the_closed_range() {
        let ceiling = InbreedingF::try_new(1.0).unwrap_err().to_string();
        let out_of_range = InbreedingF::try_new(1.5).unwrap_err().to_string();
        for message in [&ceiling, &out_of_range] {
            assert!(
                !message.contains("[0, 1]"),
                "an inbreeding rejection must not name the closed range: {message}"
            );
        }
        assert!(
            ceiling.contains("heterozygote"),
            "the ceiling message must say what F = 1 costs: {ceiling}"
        );
        assert!(
            !out_of_range.contains("heterozygote"),
            "the out-of-range message is for a typo and needs no model talk: {out_of_range}"
        );
    }

    /// Each rate names **its own quantity** when it rejects, so a message cannot
    /// send a reader to the wrong fit. `GenotypeFrequency`, `InbreedingF` and
    /// `ExpectedHeterozygosity` have their own variants; `ErrorRate` shares
    /// `DomainError::ErrorRate` with the two emission models, which is deliberate
    /// reuse — all three mean "a per-base error rate that is not a probability".
    ///
    /// **Both directions, for all four.** Each range check is two comparisons, and
    /// a test that only ever crosses one of them leaves the other free to be
    /// widened: `InbreedingF` accepting `1.5` is the live hazard, since a user
    /// types that one at a shell.
    ///
    /// `InbreedingF` has a third rejection the other two do not: exactly `1`, its
    /// excluded ceiling, which carries its own variant
    /// ([`DomainError::InbreedingFAtCeiling`]) because it is a fraction and the
    /// other message would be false of it.
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
        assert_eq!(
            InbreedingF::try_new(1.0),
            Err(DomainError::InbreedingFAtCeiling(1.0))
        );
        assert_eq!(
            ExpectedHeterozygosity::try_new(-0.5),
            Err(DomainError::ExpectedHeterozygosity(-0.5))
        );
        assert_eq!(
            ExpectedHeterozygosity::try_new(1.5),
            Err(DomainError::ExpectedHeterozygosity(1.5))
        );
        assert_eq!(
            ExpectedAlternativeFrequency::try_new(-0.5),
            Err(DomainError::ExpectedAlternativeFrequency(-0.5))
        );
        assert_eq!(
            ExpectedAlternativeFrequency::try_new(1.5),
            Err(DomainError::ExpectedAlternativeFrequency(1.5))
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
            assert!(
                matches!(ErrorRate::try_new(bad), Err(DomainError::ErrorRate(_))),
                "error rate {bad}"
            );
            assert!(
                matches!(
                    GenotypeFrequency::try_new(bad),
                    Err(DomainError::GenotypeFrequency(_))
                ),
                "genotype frequency {bad}"
            );
            assert!(
                matches!(InbreedingF::try_new(bad), Err(DomainError::InbreedingF(_))),
                "inbreeding {bad}"
            );
            assert!(
                matches!(
                    ExpectedHeterozygosity::try_new(bad),
                    Err(DomainError::ExpectedHeterozygosity(_))
                ),
                "expected heterozygosity {bad}"
            );
            assert!(
                matches!(
                    ExpectedAlternativeFrequency::try_new(bad),
                    Err(DomainError::ExpectedAlternativeFrequency(_))
                ),
                "expected alternative-allele frequency {bad}"
            );
        }
    }

    /// The species-range fallback is **one difference per thousand bases**, and the
    /// value is what is pinned.
    ///
    /// **Asserting only that it round-trips through `try_new` would pin nothing**: the
    /// constant would sit on both sides of the comparison, so every value in `[0, 1]`
    /// passes — including `0.1`, the slip that reads `1e-3` as a percentage, and `1.0`,
    /// the slip that reads it per kilobase. Both were run as mutations against this
    /// suite and both survived it. The first assertion below is what kills them.
    ///
    /// The second is the one the associated const needs: `SPECIES_FALLBACK` is built as
    /// `Self(1e-3)` and so does not pass through the range check every other value of
    /// this type does. It has to agree with the constructor, or the type has one value
    /// its own predicate never saw.
    #[test]
    fn the_species_fallback_is_one_difference_per_thousand_bases() {
        assert_eq!(ExpectedHeterozygosity::SPECIES_FALLBACK.get(), 1e-3);
        assert_eq!(
            ExpectedHeterozygosity::try_new(1e-3),
            Ok(ExpectedHeterozygosity::SPECIES_FALLBACK)
        );
    }

    /// The rejection **message** names the diversity and not a neighbouring fit, which
    /// is the whole reason the variant is separate from `GenotypeFrequency` — ng carries
    /// a heterozygosity under both types. Asserting the variant, as the tests above do,
    /// leaves the rendered text free: rewording the `#[error]` attribute to say "genotype
    /// frequency" was run as a mutation and survived the whole suite.
    #[test]
    fn expected_heterozygosity_rejection_names_its_own_quantity() {
        assert_eq!(
            ExpectedHeterozygosity::try_new(1.5)
                .unwrap_err()
                .to_string(),
            "expected heterozygosity 1.5 is not a finite probability in [0, 1]"
        );
    }

    /// `-0.0` is a probability, constructs, and comes back with the sign bit it went in
    /// with — this type is a transparent wrapper, unlike [`Phred`], which normalises the
    /// sign of zero on purpose so a quality of zero has one spelling.
    ///
    /// **Bits, not `assert_eq!`**: `-0.0 == 0.0` is true, so an accessor that quietly
    /// took an absolute value would pass an equality assertion. It is reachable rather
    /// than academic — the fitted density's segregating mass is a `.max(0.0)`
    /// (`parameter_estimation::joint`), and a product with a negative zero keeps the
    /// sign.
    #[test]
    fn expected_heterozygosity_carries_negative_zero_verbatim() {
        let zero = ExpectedHeterozygosity::try_new(-0.0).unwrap().get();
        assert_eq!(zero.to_bits(), (-0.0f64).to_bits());
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
        /// Over the whole `f64` line: each of the four rates accepts a value
        /// **exactly when** it is a finite number in its own range, and an accepted
        /// value comes
        /// back bit for bit — no clamp, no round.
        ///
        /// **Two ranges, not one.** `ErrorRate`, `GenotypeFrequency` and
        /// `ExpectedHeterozygosity` are closed `[0, 1]`; `InbreedingF` is half-open
        /// `[0, 1)`, because `F = 1` makes every heterozygote impossible
        /// (`spec/calling_priors.md` §7). Sharing one expectation across all four is
        /// what would let the ceiling be dropped again without a test noticing.
        ///
        /// This is what the point assertions above cannot do. Each range check is
        /// two comparisons, and a widened bound on either side is a value reaching
        /// step 8's genotype prior instead of an error — a wrong genotype rather
        /// than a failure. The dense `-2.0..3.0` arm is load-bearing: sampling
        /// `f64::ANY` alone never lands in `(1, 2]`, so it does not notice an
        /// `InbreedingF` whose upper bound has moved to 2.
        ///
        /// **The last two arms are the boundary itself, and they are there because
        /// neither of the first two reaches it.** Measured on the two-arm generator:
        /// a million draws produced `1.0` exactly zero times, and came no closer
        /// than 2.6 in a million. `f64::ANY` fills a `u64` and masks it, so one
        /// named bit pattern arrives about once in 2⁶⁴ draws; the dense arm halves
        /// its interval down to adjacent floats, which is uniform in real measure
        /// and so lands on any one of them about once in 10¹⁷. Without these arms
        /// the half-open expectation above would be asserted at every value except
        /// the one it was written for — green, and blind. `1.0` must be rejected
        /// and the `f64` immediately below it accepted, and the pair is what makes
        /// a moved bound a failure rather than a silent pass.
        #[test]
        fn the_constrained_rates_accept_exactly_the_probabilities_and_round_trip(
            x in proptest::prop_oneof![
                20 => proptest::num::f64::ANY,
                20 => -2.0f64..3.0f64,
                1 => proptest::strategy::Just(1.0f64),
                1 => proptest::strategy::Just(f64::from_bits(1.0f64.to_bits() - 1)),
            ]
        ) {
            let is_probability = x.is_finite() && (0.0..=1.0).contains(&x);
            let is_inbreeding_coefficient = x.is_finite() && (0.0..1.0).contains(&x);
            for (accepted, expected) in [
                (ErrorRate::try_new(x).map(ErrorRate::get), is_probability),
                (
                    GenotypeFrequency::try_new(x).map(GenotypeFrequency::get),
                    is_probability,
                ),
                (
                    InbreedingF::try_new(x).map(InbreedingF::get),
                    is_inbreeding_coefficient,
                ),
                (
                    ExpectedHeterozygosity::try_new(x).map(ExpectedHeterozygosity::get),
                    is_probability,
                ),
            ] {
                proptest::prop_assert_eq!(accepted.is_ok(), expected, "x = {}", x);
                if let Ok(value) = accepted {
                    proptest::prop_assert_eq!(value.to_bits(), x.to_bits(), "x = {}", x);
                }
            }
        }

        /// Over the whole `f32` line: a [`Phred`] is accepted **exactly when**
        /// the value is finite and at or above zero, and an accepted value comes
        /// back **bit for bit** — the one normalisation being `-0.0`, which
        /// becomes `+0.0` so a quality of zero has a single spelling.
        ///
        /// The bit-for-bit half is what the point assertions cannot do:
        /// `assert_eq!` on floats cannot separate `+0.0` from `-0.0`, so a sign
        /// flip on zero, a clamp or a round hides from them. The dense
        /// `-1.0..1.0` arm is load-bearing for the same reason the rates' arm is
        /// — sampling `f32::ANY` alone essentially never lands next to the one
        /// boundary this type exists to defend.
        #[test]
        fn phred_accepts_exactly_the_finite_non_negative_values_and_round_trips(
            q in proptest::prop_oneof![proptest::num::f32::ANY, -1.0f32..1.0f32]
        ) {
            let accepted = Phred::try_new(q).map(Phred::get);
            proptest::prop_assert_eq!(accepted.is_ok(), q.is_finite() && q >= 0.0, "q = {}", q);
            if let Ok(value) = accepted {
                let expected = if q == 0.0 { 0.0f32 } else { q };
                proptest::prop_assert_eq!(value.to_bits(), expected.to_bits(), "q = {}", q);
            }
        }

        /// [`Genotype::new`] establishes the one spelling the type rests on: the
        /// alleles come back **sorted**, they are the **same multiset** that went
        /// in, and any other order of the same alleles is the same value.
        ///
        /// Over the whole `u16` id range and up to eight copies, because every
        /// point test uses ids `0..=3` and at most four copies — and two wrong
        /// sorts reproduce those fixtures exactly. One keys on a narrowed width
        /// (`|a| a.0 as u8`) and misorders any id at or above 256; one guards the
        /// sort with a cheap `first() > last()` test for "already sorted" and
        /// leaves interior disorder in place from three copies up, which is a
        /// triploid or tetraploid heterozygote counted twice in a cohort. The
        /// dense `0..4` arm is load-bearing for the reason the rates' arm is:
        /// sampling the full `u16` alone essentially never repeats an id, and a
        /// repeated id is what homozygous means.
        #[test]
        fn a_genotype_sorts_its_alleles_and_keeps_the_multiset(
            ids in proptest::collection::vec(
                proptest::prop_oneof![proptest::num::u16::ANY, 0u16..4],
                1..=8usize,
            ),
            rotation in 0usize..8,
        ) {
            let alleles: Vec<AlleleId> = ids.iter().copied().map(AlleleId).collect();
            let genotype = Genotype::new(alleles.clone());

            proptest::prop_assert!(
                genotype.alleles().windows(2).all(|pair| pair[0] <= pair[1]),
                "not sorted: {:?}",
                genotype.alleles()
            );

            let mut same_multiset = alleles.clone();
            same_multiset.sort_unstable();
            proptest::prop_assert_eq!(
                genotype.alleles(),
                &same_multiset[..],
                "an allele was dropped, added or altered"
            );

            let mut respelled = alleles.clone();
            let copies = respelled.len();
            respelled.rotate_left(rotation % copies);
            proptest::prop_assert_eq!(
                &genotype,
                &Genotype::new(respelled),
                "two spellings of one genotype must be one value"
            );
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

    /// A batch renders as the bare index, so a message can supply its own word for it —
    /// "batch {batch}", not "batch batch 2". Untested until mutation testing gave `Display` a
    /// prefix and the whole suite stayed green (C2's review).
    #[test]
    fn a_batch_renders_as_the_index_alone() {
        assert_eq!(BatchId(2).to_string(), "2");
        assert_eq!(BatchId::ALL_TOGETHER.to_string(), "0");
    }

    /// The default batching is batch zero, which is what makes an all-zero batching mean
    /// "every read group ran together" without anyone saying so.
    #[test]
    fn the_default_batch_is_the_first_one() {
        assert_eq!(BatchId::ALL_TOGETHER.get(), 0);
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

    /// The boundary in both directions: a quality of exactly zero is legal — it
    /// is `p = 1`, a call that cannot be wrong — and anything below it is not,
    /// because a negative Phred is a probability above one. `NaN` and `-∞` go
    /// the same way as the negative, since `quality >= 0.0` is false for both.
    /// `+∞` is rejected too but as a different event, so it has its own variant.
    ///
    /// The top of the range is deliberately open: this type refuses to cap,
    /// because where to cap a `GQ` is the consumer's decision (see the type's
    /// doc comment), so a ceiling appearing inside `try_new` later must break a
    /// test.
    #[test]
    fn phred_accepts_zero_and_rejects_everything_below_it() {
        assert_eq!(Phred::try_new(0.0).unwrap().get(), 0.0);
        assert_eq!(Phred::try_new(30.0).unwrap().get(), 30.0);
        assert_eq!(Phred::try_new(f32::MAX).unwrap().get(), f32::MAX);
        assert!(matches!(
            Phred::try_new(-f32::EPSILON),
            Err(DomainError::Phred(q)) if q < 0.0
        ));
        assert_eq!(Phred::try_new(-1.0), Err(DomainError::Phred(-1.0)));
        // Comparison with `NaN` is never true, so the `>= 0.0` test rejects it
        // — and `DomainError`'s `PartialEq` is IEEE equality on the payload, so
        // this must be a `matches!` and not an `assert_eq!`.
        assert!(matches!(
            Phred::try_new(f32::NAN),
            Err(DomainError::Phred(q)) if q.is_nan()
        ));
        assert_eq!(
            Phred::try_new(f32::INFINITY),
            Err(DomainError::PhredInfinite)
        );
        assert!(matches!(
            Phred::try_new(f32::NEG_INFINITY),
            Err(DomainError::Phred(q)) if q == f32::NEG_INFINITY
        ));
    }

    /// A quality of zero has one spelling, whichever door it came in by.
    ///
    /// `-PHRED_PER_NAT * 0.0` is `-0.0` under IEEE, and `-0.0` passes
    /// `>= 0.0` — it *is* zero — so nothing in the range check notices. It
    /// matters because this type's purpose is VCF's `QUAL` and `GQ`, where a
    /// certainty printing as `-0` is not a number those columns should carry.
    /// `assert_eq!(.., 0.0)` cannot see this: `-0.0 == 0.0` is true.
    #[test]
    fn phred_zero_is_positive_zero_whichever_constructor_made_it() {
        let from_certainty = Phred::from_log_prob(LogProb(0.0))
            .expect("p = 1 is quality zero")
            .get();
        assert!(
            from_certainty.is_sign_positive(),
            "quality zero must be +0.0, not -0.0 (bits {:#x})",
            from_certainty.to_bits()
        );
        assert_eq!(
            format!("{from_certainty}"),
            "0",
            "a QUAL column cannot say -0"
        );
        assert!(Phred::try_new(-0.0).unwrap().get().is_sign_positive());
    }

    /// The conversion against numbers worked out by hand rather than by the
    /// same formula: `-10 log10(p)` is 30 at one wrong call in a thousand and
    /// 20 at one in a hundred, which is what those Phreds mean, and
    /// `-10 log10(0.5)` is `10 log10 2` = 3.0103.
    ///
    /// **The two decade pairs are the ones that discriminate**, because the
    /// `1e-4` tolerance is absolute and these values are an order of magnitude
    /// apart: at 30 it admits a relative error of 3.3e-6 in the scale factor, at
    /// 3.0103 ten times that. What sits outside it at 30: an implementation
    /// using `log2` instead of `log10` returns 99.657845 where this asserts 30,
    /// one that dropped the factor of ten returns exactly 3, and one using
    /// `baq`'s htslib-parity `4.343` returns 30.000381. The `0.5` pair earns its
    /// place only as the one expected value that is not a round number.
    #[test]
    fn phred_from_log_prob_matches_the_hand_computed_scale() {
        let phred_of = |p: f64| Phred::from_log_prob(LogProb(p.ln())).unwrap().get();
        assert!((phred_of(0.001) - 30.0).abs() < 1e-4, "{}", phred_of(0.001));
        assert!((phred_of(0.01) - 20.0).abs() < 1e-4, "{}", phred_of(0.01));
        assert!((phred_of(0.5) - 3.010_3).abs() < 1e-4, "{}", phred_of(0.5));
    }

    /// What the doc comment on [`Phred::from_log_prob`] promises about widths:
    /// the scaling happens at `f64` and narrows once, at the end.
    ///
    /// Pinned at a quality of 3000 because that is where the ordering shows.
    /// Narrowing `ln p` first and multiplying in `f32` returns 2999.9998 here,
    /// one `f32` step low, while the three everyday qualities of the test above
    /// (30, 20, 3.0103) come out identical either way — so without this the
    /// invariant has no test behind it.
    #[test]
    fn phred_from_log_prob_keeps_full_f64_width_before_narrowing() {
        let quality = Phred::from_log_prob(LogProb(1e-300f64.ln())).expect("a finite quality");
        assert_eq!(
            quality.get().to_bits(),
            3000.0f32.to_bits(),
            "got {} (bits {:#x})",
            quality.get(),
            quality.get().to_bits()
        );
    }

    /// The three things the scale cannot hold, told apart by the error rather
    /// than by the caller inspecting a float. `ln(1) = 0` is the one that must
    /// NOT be an error: a call that cannot be wrong is quality zero, the low end
    /// of the scale.
    #[test]
    fn phred_from_log_prob_rejects_what_the_scale_cannot_hold() {
        // p = 1 — certainty. Quality zero, and legal.
        assert_eq!(Phred::from_log_prob(LogProb(0.0)).unwrap().get(), 0.0);
        // p = 0 — the score of an impossible read line-up, which `LogProb`
        // carries deliberately. Its quality is infinite, which Phred does not
        // hold; capping it is the consumer's decision, not this type's, and the
        // variant is what tells that consumer this was not a broken sum.
        assert_eq!(
            Phred::from_log_prob(LogProb(f64::NEG_INFINITY)),
            Err(DomainError::PhredInfinite)
        );
        // A log probability finite in `f64` whose scaled quality saturates
        // `f32`. Unreachable from real data — it needs `ln p` below about
        // -7.8e37 — and pinned so the `is_finite` guard is not simplified away.
        assert_eq!(
            Phred::from_log_prob(LogProb(-1e300)),
            Err(DomainError::PhredInfinite)
        );
        // A positive logarithm is a probability above one — the caller's
        // arithmetic is wrong, and the negative quality says so.
        assert!(matches!(
            Phred::from_log_prob(LogProb(0.1)),
            Err(DomainError::Phred(q)) if q < 0.0
        ));
        // `LogProb`'s field is public and unconstrained, so a caller's
        // `0.0 / 0.0` — an unnormalised posterior, say — arrives here as a
        // `NaN`. It must not reach a `QUAL` column as a silently wrong record.
        assert!(matches!(
            Phred::from_log_prob(LogProb(f64::NAN)),
            Err(DomainError::Phred(q)) if q.is_nan()
        ));
    }

    /// The property the whole design of this type rests on: a genotype is a
    /// **multiset**, so the order the alleles were handed over in cannot change
    /// what it equals.
    ///
    /// It matters because the calling loop mints one genotype per sample from a
    /// table row and the caller downstream groups, compares and hashes them. If
    /// `[1, 0]` and `[0, 1]` were different values, one heterozygote in a cohort
    /// would count as two, and every derived `Hash` and `Ord` would disagree
    /// with the `PartialEq` a reader assumes.
    #[test]
    fn a_genotypes_alleles_are_a_multiset_not_a_sequence() {
        use std::cmp::Ordering;

        let reference_first = Genotype::new(vec![AlleleId(0), AlleleId(1)]);
        let alternate_first = Genotype::new(vec![AlleleId(1), AlleleId(0)]);
        assert_eq!(reference_first, alternate_first);
        // Sorted, so `alleles()` reads the same either way — which is what makes
        // the derived `Ord` and `Hash` agree with that equality.
        assert_eq!(reference_first.alleles(), alternate_first.alleles());
        assert_eq!(reference_first.alleles(), [AlleleId(0), AlleleId(1)]);

        // Both follow from the assertion above, since `Genotype` has one field
        // and `PartialEq`, `Hash` and `Ord` are all derived from it. They are
        // here to fail if that ever stops being true — a second field, or a
        // hand-written impl of any of the three.
        let mut distinct_genotypes = std::collections::HashSet::new();
        distinct_genotypes.insert(reference_first.clone());
        distinct_genotypes.insert(alternate_first.clone());
        assert_eq!(
            distinct_genotypes.len(),
            1,
            "one genotype, however it was spelled"
        );
        assert_eq!(reference_first.cmp(&alternate_first), Ordering::Equal);
    }

    /// A multiset, not a set: two copies of one allele is what homozygous means,
    /// so the repeat must survive construction. A `sort` that deduplicated, or a
    /// `HashSet` reached for because "order does not matter", would turn every
    /// homozygote into a haploid call and nothing about the type would object.
    #[test]
    fn a_genotype_keeps_repeated_alleles() {
        let homozygous_alternate = Genotype::new(vec![AlleleId(1), AlleleId(1)]);
        assert_eq!(
            homozygous_alternate.alleles(),
            [AlleleId(1), AlleleId(1)],
            "both copies are carried"
        );
        assert_ne!(
            homozygous_alternate,
            Genotype::new(vec![AlleleId(1)]),
            "a diploid homozygote is not a haploid call"
        );
    }

    /// How many alleles a genotype holds is the ploidy at that region. The
    /// haploid and tetraploid cases pin that `new` neither pads nor truncates,
    /// and that sorting works past the two entries every other point test uses.
    ///
    /// **"No ceiling on ploidy" is not what these two points show, and no point
    /// fixture could show it** — that is a claim about every length, and unlike
    /// [`Ploidy`], whose domain is a finite `u8` and so can be enumerated, a
    /// genotype's length domain is not. The property test reaches eight copies;
    /// past that the guarantee rests on `Box<[AlleleId]>` having no length limit
    /// of its own, not on a test.
    ///
    /// The tetraploid fixture **starts and ends in order** on purpose. Written
    /// largest-first and smallest-last it would trip any "is this already
    /// reversed?" fast path into sorting anyway, so it would pass against a
    /// `new` that skipped the sort whenever the first entry was not above the
    /// last — which leaves interior disorder in place from three copies up.
    #[test]
    fn a_genotype_holds_one_allele_per_genome_copy() {
        let haploid = Genotype::new(vec![AlleleId(2)]);
        assert_eq!(haploid.alleles().len(), 1);

        let tetraploid = Genotype::new(vec![AlleleId(0), AlleleId(3), AlleleId(2), AlleleId(0)]);
        assert_eq!(tetraploid.alleles().len(), 4);
        assert_eq!(
            tetraploid.alleles(),
            [AlleleId(0), AlleleId(0), AlleleId(2), AlleleId(3)]
        );
    }

    /// An empty multiset is the one construction that is not a genotype at any
    /// ploidy: [`Ploidy`] refuses zero copies because a genome with none is not
    /// a genome, and `alleles().len()` **is** that ploidy. Today nothing can
    /// reach it — the only minter expands a genotype-table row, and a row holds
    /// exactly `ploidy` copies — which is the reason to pin it now rather than
    /// after a row builder learns to emit a zero-length row.
    #[test]
    #[should_panic(expected = "one allele per genome copy")]
    fn a_genotype_cannot_be_built_from_no_alleles_at_all() {
        let _ = Genotype::new(vec![]);
    }

    /// Re-minting a genotype from what `alleles()` handed out gives the same
    /// value back. `new` canonicalises, and a canonical form has to be a fixed
    /// point: a rule that sorted and then rotated would still make the two
    /// spellings in `a_genotypes_alleles_are_a_multiset_not_a_sequence` agree
    /// with each other, while changing a genotype every time it was rebuilt —
    /// which is what a caller does when it widens a diploid call, or
    /// reconstructs a genotype read back from a VCF.
    #[test]
    fn a_genotype_new_is_idempotent_on_its_own_alleles() {
        for spelling in [
            vec![AlleleId(1), AlleleId(0)],
            vec![AlleleId(1), AlleleId(1)],
            vec![AlleleId(3), AlleleId(0), AlleleId(2), AlleleId(0)],
        ] {
            let once = Genotype::new(spelling);
            let twice = Genotype::new(once.alleles().to_vec());
            assert_eq!(once, twice, "canonical form must be a fixed point");
            assert_eq!(once.alleles(), twice.alleles());
        }
    }

    // -----------------------------------------------------------------
    // SummedLogError
    // -----------------------------------------------------------------

    /// The step is what a psp's header records so a reader can interpret the integer, so the
    /// two directions have to be each other's inverse on the grid.
    #[test]
    fn a_summed_log_error_round_trips_through_its_own_step() {
        for steps in [0i64, 1, -1, -61_030, 4_096, -8_388_608] {
            let value = SummedLogError::from_steps(steps);
            assert_eq!(value.steps(), steps);
            assert_eq!(SummedLogError::from_nats(value.nats()), value);
        }
    }

    /// **The error a step introduces is the whole risk of the type**, since this term goes
    /// straight into a genotype likelihood. Half a step is 1.22 × 10⁻⁴ natural logs, which the
    /// spec prices at 0.024 % of the likelihood term.
    #[test]
    fn rounding_never_moves_a_value_by_more_than_half_a_step() {
        let half_a_step = 0.5 / SummedLogError::STEPS_PER_NAT as f64;
        // Values chosen to land at, just under and just over a step boundary, and at the
        // depths the caller actually sees: −3 at three reads a position, −3360 at three
        // hundred (the locus D3's real-data run found).
        for nats in [
            0.0,
            -1e-9,
            -0.5,
            -2.999_999_999_999_999,
            -3.0,
            -14.9,
            -3_360.392_684_715,
        ] {
            let rounded = SummedLogError::from_nats(nats).nats();
            assert!(
                (rounded - nats).abs() <= half_a_step,
                "{nats} rounded to {rounded}, which is more than half a step away"
            );
        }
    }

    /// **The property the type exists for.** Two routes that add the same reads' error in
    /// different orders must reach the same value — which `f64` addition does not guarantee and
    /// integer addition does.
    #[test]
    fn adding_is_exact_and_does_not_care_about_order() {
        let terms: Vec<SummedLogError> = [-0.1, -7.6, -0.000_3, -12.5, -3.0]
            .iter()
            .map(|&nats| SummedLogError::from_nats(nats))
            .collect();

        let forwards: SummedLogError = terms.iter().copied().sum();
        let backwards: SummedLogError = terms.iter().rev().copied().sum();
        assert_eq!(forwards, backwards);

        // The same addends in `f64` do not have that property, which is why this type is not
        // an `f64`. One large term and two that fall below its last bit: added large-first the
        // small ones vanish, added small-first they survive to move the total. If this
        // assertion ever fails, the fixture stopped exercising the point.
        const UNEVEN: [f64; 3] = [-1.0, -1e-16, -1e-16];
        let sum_forwards: f64 = UNEVEN.iter().sum();
        let sum_backwards: f64 = UNEVEN.iter().rev().sum();
        assert_ne!(
            sum_forwards, sum_backwards,
            "the fixture must contain addends that `f64` adds differently in each order"
        );
    }

    /// A not-a-number error sum once came back through `f64::max` as the most confident read
    /// the model can express — a confident wrong answer with nothing failing. The type has no
    /// such value, and the conversion says what it does instead of leaving it to `as`.
    #[test]
    fn a_value_no_read_can_produce_becomes_a_stated_one_rather_than_a_silent_one() {
        assert_eq!(SummedLogError::from_nats(f64::NAN), SummedLogError::NONE);
        assert_eq!(
            SummedLogError::from_nats(f64::NEG_INFINITY).steps(),
            i64::MIN
        );
        assert_eq!(SummedLogError::from_nats(f64::INFINITY).steps(), i64::MAX);
        // And saturating rather than wrapping, so an extreme value stays extreme.
        let most = SummedLogError::from_steps(i64::MAX);
        assert_eq!((most + most).steps(), i64::MAX);
    }

    /// The unit is in the rendering, because a bare number here reads as a probability, a
    /// Phred score or a step count depending on who is looking.
    #[test]
    fn a_summed_log_error_renders_with_its_unit() {
        assert_eq!(SummedLogError::from_nats(-3.0).to_string(), "-3 ln");
        assert_eq!(SummedLogError::NONE.to_string(), "0 ln");
    }
}
