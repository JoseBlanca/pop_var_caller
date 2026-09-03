//! The STR row — one sample's genotype likelihoods at one repeat tract.
//!
//! *How well does each candidate genotype explain what this sample's reads showed at this
//! tract?* The answer is one number per genotype, and this module computes the whole row.
//!
//! Everything about **how a single read is scored against a single allele** lives behind the
//! emission seam ([`super::ssr_emission`]); everything about **how those scores are combined
//! into a genotype's likelihood** lives here. That split is what lets the model comparison
//! behind spec §4.1 swap one model for another without touching the arithmetic around it.
//!
//! # The formula, in words
//!
//! **A read arrived one of three ways**: copied from one of this individual's own copies of the
//! tract, from something no candidate explains at all, or from **another individual's DNA in
//! the same library**. Writing that for a genotype `g` whose copy counts are `k_a` over a
//! ploidy `P` (spec §2.1, §4.5, §4.5.1):
//!
//! ```text
//! log Lg(g)  =  Σ_o  n_o · log[ (1 − λ − c) · Σ_a (k_a / P) · Lr(o | a)
//!                             +      λ      · U
//!                             +      c      · seed(o) ]
//! ```
//!
//! - `o` runs over the sample's observations at this tract and `n_o` is how many reads showed
//!   each one;
//! - `k_a / P` is the chance a read was copied from a copy carrying allele `a`;
//! - `Lr(o | a)` is the emission — [`SsrEmissionModel::emission`] for a read that spanned the
//!   tract, [`SsrEmissionModel::censored_emission`] for one that ran off its own end;
//! - `λ` is the chance the read explains nothing, and `U` is what such a read shows: uniform
//!   over the tract lengths the candidates can reach;
//! - `c` is the share of this **read group's** reads that came from somebody else, and
//!   `seed(o)` is how common the length that observation showed is at this locus.
//!
//! **The second term is what stops one strange read from ruling out every genotype.** A read
//! from a paralogous tract, a chimera, a somatic length in a long tract — its emission is zero
//! under every candidate, so the bracket collapses to `λ · U` for every genotype alike and the
//! term contributes the same number to every entry of the row.
//!
//! **The third is a different thing and folding it into the second would be wrong.** A junk
//! read can show anything, so its distribution is flat and its term cancels; a contaminating
//! read shows a *plausible* length, so its distribution is peaked and its term does not.
//! Adding `c` to `λ` and keeping one flat distribution would give a contaminant's reads the
//! junk treatment — which is what the model does today, and what spec §4.5.1 exists to stop.
//! With `contamination: None` the row is the two-term form.
//!
//! # What a row costs, and why that is the design
//!
//! **`observations × candidates` emission calls, not `observations × genotypes`.** Each
//! emission is computed once into [`SsrRowScratch`]'s cache and read by every genotype that
//! carries the candidate. At six candidates and a diploid that is 6 calls an observation
//! rather than 21 — and spec §8 calls it the design rather than an optimisation, which is why
//! the cache is a field of the caller's scratch rather than a local here.
//!
//! # Where the two distributions come from
//!
//! **The lengths `U` is spread over are the locus's, not the cohort's**
//! ([`super::ssr_emission::fill_reachable_lengths`]). Production spreads its outlier weight
//! over the number of distinct sequences the whole cohort showed, which made a sample's junk
//! floor ten times lower in a 63-accession panel than alone (spec §4.5); this asks the
//! candidates and the two slip cutoffs and nothing else.
//!
//! **The seed comes from the genotype prior's starting shape**, converted to a distribution
//! over those same lengths by the calling loop — which is the only place holding both the
//! candidate table and this support. Spec §4.5.1 says why it is the prior's shape and not the
//! cohort's fitted frequencies: the second would move on every pass, and contamination is
//! frozen.

use super::ContaminationView;
use super::ssr_emission::{SsrCandidate, SsrEmissionModel, SsrScoringContext};
use super::{SsrRowScratch, SsrSampleEvidence};
use crate::ng::calling::GenotypeTableView;
use crate::ng::locus_generation::{ReadWitness, SequenceObservation};
use crate::ng::parameter_estimation::Provenance;
use crate::ng::types::{LogProb, ReadGroupId};

/// How often a read at a repeat tract came from somewhere other than this individual's copies
/// of it — **0.05, chosen by a sweep against genotype accuracy, and it is not a measurement of
/// that share.**
///
/// # What it does, which is not what its name says
///
/// The junk term spreads this weight evenly over every tract length the model can reach — about
/// twenty-two of them at a homopolymer — so every read's emission has a floor of roughly
/// `weight / lengths` under it, whatever genotype is being scored. **A floor on the emission is
/// a cap on how much one read can pull a genotype**: past it, a read being more surprising buys
/// no more evidence. freebayes does that job with a read-dependence factor and GATK with a
/// Phred-45 cap; this is the only thing in ng doing it at a tract.
///
/// **Where the floor sits against the stutter distribution is the whole effect.** At 0.01 the
/// floor is 4.6 x 10^-4, five times *below* the chance of a read slipping two whole repeats
/// (2.4 x 10^-3), so such a read scores as real evidence for a second allele. At 0.05 the floor
/// is 2.3 x 10^-3 — level with that slip — so a two-repeat slip product carries almost no
/// evidence either way, which is the intended behaviour and the change that moves the numbers
/// below.
///
/// # Why 0.05, and what that number is worth
///
/// **Measured**, on GIAB's HG002 tandem-repeat benchmark at 30x, 20,204 typed tracts, one full
/// run a setting scored against the assembly-based truth
/// (`doc/devel/reports/ng_tract_genotype_improvement_2026-09-02.md` §5.2):
///
/// | weight | homopolymer | period 2+ | heterozygote called for a homozygous truth |
/// |---|---|---|---|
/// | 0.01, the inherited value | 0.8851 | 0.9037 | 88 |
/// | **0.05** | **0.8881** | **0.9059** | 77 |
/// | 0.10 | 0.8892 | 0.9059 | 73 |
/// | 0.20 | 0.8887 | 0.9051 | 70 |
/// | 0.30 | 0.8891 | 0.9043 | 63 |
///
/// The curve is flat from 0.05 to 0.30 at homopolymers and falls away above 0.10 at period 2 and
/// above, so what the sweep really says is **"not 0.01"**. 0.05 is the owner's choice of the
/// conservative end of that plateau (2026-09-03): it takes about three quarters of the available
/// gain while moving least far from the value the reads themselves suggest.
///
/// # What is still open, and why this is not called fitted
///
/// **The literal reading of this number disagrees with the sweep.** Read as what it is named —
/// the share of reads nothing explains — it measures 1 in 2,300 at homopolymers and 1 in 209 at
/// period 2 and above, which is both far below 0.05 and ordered the opposite way to what the
/// sweep prefers. So the constant is doing a job nobody named it for, and its warrant stays
/// `Defaulted`: a stated constant, not an estimate.
///
/// **Three things it has not been tested against**, all of them inside the range this caller is
/// committed to (`doc/devel/ng/spec/design_principles.md` §0): a second individual, a cohort, and
/// low depth. A floor under every read's emission behaves very differently at three reads a tract
/// than at thirty, and the sweep was run at 30x and 50x on one sample. **Sweeping it per motif
/// period, and on the tomato panel at three reads, is the work this constant owes.**
///
/// Production's value, which ng carried until now, is 0.01
/// ([`em.rs`](../../../../src/ssr/cohort/em.rs)).
pub const DEFAULT_OUTLIER_WEIGHT: f64 = 0.05;

/// **The outlier weight this run scored with, and whether the run was handed it or inherited
/// it.**
///
/// **Two states and no more.** Nothing fits this number, so [`Provenance::FittedHere`] and
/// [`Provenance::Borrowed`] are unreachable for it: either the run read a value out of a
/// parameters file, which is `Supplied`, or it took [`DEFAULT_OUTLIER_WEIGHT`], which is
/// `Defaulted`. **Both fields are private and the two constructors below set them together**,
/// so no caller outside this module can write the unreachable pair, and neither constructor
/// can be given half of one — [`Self::defaulted`] fixes the value as well as the warrant.
///
/// **The file can spell the other two and the reader refuses them.** `Warrant` in the
/// parameters file has all four states because most of its numbers need all four, so
/// `ParametersFile::validate` is what keeps this key to the two — including a `defaulted`
/// value that is not the compiled-in constant, which is what a person who edits the number
/// and leaves the warrant alone produces.
///
/// **Why it is a pair rather than an `f64`.** `doc/devel/ng/spec/parameters_file.md` §3.8 puts
/// this number in the file *so that a person can change it* — "marking it soft is the point of
/// writing it down" — and a run that kept only the number could not tell an edited 0.01 from
/// the compiled-in one. The file it writes afterwards would then mark a supplied value
/// `defaulted`, which is the file, the report and the score disagreeing while looking wired up.
///
/// **It is reported once for the run and never folded into a repeat tract's per-cell warrant**
/// ([`RunParameterReport::repeat_tract_outlier_weight`](crate::ng::calling::run_report::RunParameterReport::repeat_tract_outlier_weight)):
/// it is one run-wide number, so folding it in would mark *every* tract of every run as
/// resting on a defaulted parameter — or, under a supplied weight, on a supplied one, which the
/// ladder ranks only a rung above — and erase the fitted-against-borrowed distinction the
/// per-cell warrant exists to carry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RepeatTractOutlierWeight {
    value: f64,
    provenance: Provenance,
}

impl RepeatTractOutlierWeight {
    /// [`DEFAULT_OUTLIER_WEIGHT`], marked `Defaulted` — what a run that was handed no value gets.
    #[must_use]
    pub const fn defaulted() -> Self {
        Self {
            value: DEFAULT_OUTLIER_WEIGHT,
            provenance: Provenance::Defaulted,
        }
    }

    /// A value the run was handed — from a parameters file — marked `Supplied`.
    ///
    /// # Panics
    ///
    /// Unless `value` is finite and strictly inside 0 and 1, which is the range
    /// [`SsrLocusParameters`] already refuses outside of, several frames later and naming a
    /// locus rather than the file the number came from. A weight of zero says no read at a
    /// repeat tract can have come from anywhere but this individual's two copies; a weight of
    /// one says none of them came from there.
    ///
    /// **A caller reading a validated parameters file does not reach this panic**:
    /// `ParametersFile::validate` refuses this key outside the same open interval, naming it.
    /// This is the guard for a caller that did not validate.
    ///
    /// **One narrower check is still the scoring row's and cannot move here.** That row also
    /// asserts `weight + contamination fraction < 1` per read group, and the fractions are
    /// the fit's rather than the file's, so a weight of 0.9 passes both this and validation
    /// and fails at the first locus of a contaminated run.
    #[must_use]
    pub fn supplied(value: f64) -> Self {
        assert!(
            value.is_finite() && value > 0.0 && value < 1.0,
            "a repeat-tract outlier weight is a share of reads strictly inside 0 and 1, and \
             {value} was supplied; a zero says no read at a tract can have come from anywhere \
             but this individual's copies of it, and a one says none of them did"
        );
        Self {
            value,
            provenance: Provenance::Supplied,
        }
    }

    /// The number itself — what a scoring row charges the junk term.
    #[inline]
    #[must_use]
    pub fn value(self) -> f64 {
        self.value
    }

    /// Whether the run was handed it or inherited it.
    #[inline]
    #[must_use]
    pub fn provenance(self) -> Provenance {
        self.provenance
    }
}

/// The smallest share of a read the row will leave to this individual's own copies.
///
/// **`1 − λ − c` is floored positive** (spec §4.5.1). It is nowhere near this at the values in
/// play — an outlier weight of 1 in 100 against a contamination fraction of 1 to 3 in 100 — and
/// the floor exists for the parameter row that puts the two shares above one, where the
/// alternative is a negative weight and a `NaN` out of `ln`.
const MIN_SHARE_FROM_THIS_INDIVIDUAL: f64 = 1e-12;

/// The scoring parameters for every `(read group, candidate)` pair at one locus, as one
/// checked table.
///
/// **The stride has one spelling, which is the whole point of the type.** A context is built
/// per `(read group, candidate)` — a read's chance of slipping is a property of the tract it
/// was copied from, so a 6-repeat candidate and a 12-repeat one at the same locus are drawn
/// from different strata (spec §4.4) — and a bare slice with the indexing written at the call
/// site is two chances to write it differently. Getting it wrong reads a real context that
/// belongs to another candidate: a plausible number, no panic.
#[derive(Debug, Clone, Copy)]
pub struct SsrScoringContextTable<'a> {
    contexts: &'a [SsrScoringContext<'a>],
    candidates: usize,
}

impl<'a> SsrScoringContextTable<'a> {
    /// Wrap the table, checking it is rectangular.
    ///
    /// # Panics
    ///
    /// If the slice is not `read_groups × candidates` entries. A short table would otherwise
    /// surface at whichever locus first reached past its end, or never.
    #[must_use]
    pub fn new(contexts: &'a [SsrScoringContext<'a>], candidates: usize) -> Self {
        assert!(
            candidates > 0,
            "a locus is called over at least its reference allele"
        );
        assert!(
            contexts.len().is_multiple_of(candidates),
            "a context table holds one entry per (read group, candidate): {} entries is not a \
             whole number of rows of {candidates}",
            contexts.len()
        );
        Self {
            contexts,
            candidates,
        }
    }

    /// How many read groups this table covers.
    #[must_use]
    pub fn read_group_count(&self) -> usize {
        self.contexts.len() / self.candidates
    }

    /// How many candidates each read group's row covers.
    #[must_use]
    pub fn candidate_count(&self) -> usize {
        self.candidates
    }

    /// The context for one `(read group, candidate)`.
    ///
    /// # Panics
    ///
    /// **In release as well as debug**, on a read group or candidate past what the table
    /// covers — because the alternative is scoring a read against another candidate's stutter
    /// parameters and getting a number back.
    #[must_use]
    pub fn of(&self, read_group: ReadGroupId, candidate: usize) -> &SsrScoringContext<'a> {
        let group = read_group.get() as usize;
        assert!(
            candidate < self.candidates,
            "candidate {candidate} is past the {} this table covers",
            self.candidates
        );
        assert!(
            group < self.read_group_count(),
            "read group {group} is past the {} this table covers",
            self.read_group_count()
        );
        &self.contexts[group * self.candidates + candidate]
    }
}

/// **A read from another individual's DNA** — spec §4.5.1's third term, which is not the junk
/// term and must not be folded into it.
///
/// A junk read can show anything, so its distribution is flat and its term is the same number
/// under every genotype; **a contaminating read shows a *plausible* length**, so its
/// distribution is peaked and its term does not cancel. Adding `c` to `λ` and keeping one flat
/// distribution would give a contaminant's reads the junk treatment, which is what the model
/// does today and what this term exists to stop.
#[derive(Debug, Clone, Copy)]
pub struct SsrContaminationMixture<'a> {
    /// How much of each read group's DNA came from somebody else, in read-group order.
    ///
    /// **Per read group and not one number**, because that is what the pre-pass fits and what
    /// spec §4.5.1 asks for: contamination is a property of a library. (The architecture's
    /// sketch wrote a scalar here; the specification's own wording is the one followed.)
    ///
    /// **[`ContaminationView`] rather than a bare `f64`**, which is the type the pre-pass
    /// produces and the SNP/indel path already stores. It carries *whose reads the fraction was
    /// fitted from* beside the number — spec §3.6's third requirement — and a second array of
    /// the same values would drop that on the floor for no gain.
    pub fraction_of_each_read_group: &'a [ContaminationView],
    /// How common each reachable tract length is — **parallel to
    /// [`SsrLocusParameters::reachable_lengths`]**, entry for entry, and summing to one.
    ///
    /// # This resolves an open question, and the resolution is worth reading
    ///
    /// The genotype prior builds its seed shape with one entry **per candidate**
    /// (`fill_seed_share_per_candidate`), and `c · seed(o)` asks for a probability per observed
    /// **length**. `arch/read_likelihoods.md` §4.1 and `arch/calling_priors.md` §5 both record
    /// the gap as open, with three cases it turns on. Keying this to the locus's reachable
    /// lengths answers all three:
    ///
    /// - **Two candidates spelling one length** — one entry, so whoever builds this sums the
    ///   two rungs rather than letting each take the full share.
    /// - **A read at a length no candidate reaches** — no entry, so the seed gives it nothing
    ///   and it falls to the outlier floor. Spec §4.5.1 says that is right rather than a loss:
    ///   *"its reads simply are not covered by the seed and fall to the outlier floor instead —
    ///   which is where they go today, so nothing is lost."*
    /// - **A read that ran out inside the tract**, which has no length at all — it gets the
    ///   seed's mass at or above what it witnessed, the same lower-bound reading §5.2 gives the
    ///   emission.
    ///
    /// **Converting the prior's per-candidate shape into this is the calling loop's job**, not
    /// this module's: the loop is what assembles a locus, and it is the only place that holds
    /// both the candidate table and this support.
    ///
    /// # Two things the resolution above does not settle
    ///
    /// - **Most entries will be zero, and that is the shape rather than a defect.** The
    ///   reachable support is far wider than the candidate set — it is every length the slip
    ///   cutoffs admit from any candidate — so a dinucleotide locus with five candidates has
    ///   about 39 reachable lengths and mass at five of them. A contaminating read at a length
    ///   no candidate carries falls to the outlier floor, which is where spec §4.5.1 puts it.
    /// - **How two candidates spelling one length should share that length's mass is still the
    ///   specification's open question 3**, not this type's. Keying to lengths is what makes
    ///   the question answerable in one place instead of being decided twice by accident;
    ///   `arch/calling_priors.md` §5 carries the measurement it turns on.
    pub contaminant_length_frequencies: &'a [f64],
}

impl SsrContaminationMixture<'_> {
    /// How much of this read group came from somebody else.
    ///
    /// # Panics
    ///
    /// On a read group past the table, in release as well as debug — the alternative is
    /// charging a read the wrong library's contamination and getting a number back.
    #[must_use]
    pub fn fraction_of(&self, read_group: ReadGroupId) -> f64 {
        let group = read_group.get() as usize;
        assert!(
            group < self.fraction_of_each_read_group.len(),
            "read group {group} is past the {} this contamination table covers",
            self.fraction_of_each_read_group.len()
        );
        self.fraction_of_each_read_group[group].fraction
    }

    /// **How likely a contaminating read is to show what this observation showed.**
    ///
    /// A read that spanned the tract showed a length, and the seed says how common that length
    /// is. A read that ran out shows a lower bound, so it gets the seed's mass at or above what
    /// it witnessed — the same reading spec §5.2 gives the emission, and for the same reason: a
    /// truncated read has not shown a short allele.
    ///
    /// A length no candidate reaches is not in the support and scores zero here, which sends
    /// the read to the outlier floor.
    fn contaminant_frequency_of(
        &self,
        observation: &SequenceObservation,
        reachable_lengths: &[u32],
    ) -> f64 {
        self.contaminant_length_frequencies
            [lengths_the_observation_allows(observation, reachable_lengths)]
        .iter()
        .sum()
    }
}

/// What the **locus** contributes to a row — the same for every sample called at this tract,
/// and none of it derivable from the genotype table.
///
/// **Grouped rather than passed loose**, and not only to shorten the signature: these five
/// travel together and are built once per locus by the calling loop. A caller that has one of
/// them has all of them.
#[derive(Debug, Clone, Copy)]
pub struct SsrLocusParameters<'a> {
    /// The candidates, in genotype-table allele order.
    ///
    /// **Built [`SsrCandidate`]s rather than bases, because a candidate's repeat count is not
    /// derivable from its bases**: an interrupted tract's byte length divided by the period is
    /// not how many repeats it holds. The locus generator has already measured it, and
    /// re-measuring here would be the duplication spec §7 puts on the alignment module's side
    /// of the boundary.
    pub candidates: &'a [SsrCandidate<'a>],
    /// The stutter and substitution parameters, per `(read group, candidate)` — never hoisted
    /// out of the candidate loop (spec §4.4).
    pub contexts: SsrScoringContextTable<'a>,
    /// How often a read came from somewhere other than this individual's copies of the tract.
    ///
    /// **The run's own value**, which is
    /// [`FrozenParameters::repeat_tract_outlier_weight`](crate::ng::calling::FrozenParameters::repeat_tract_outlier_weight)
    /// stripped of its warrant — [`DEFAULT_OUTLIER_WEIGHT`] where the run was handed none, and
    /// whatever a parameters file supplied where it was.
    pub outlier_weight: f64,
    /// The tract lengths the outlier weight is spread over — **a property of the candidate set
    /// and the two cutoffs, with no cohort in it** (spec §4.5), built by
    /// [`fill_reachable_lengths`].
    ///
    /// It is the support the contamination seed is keyed to as well, so the two terms of the
    /// mixture that are not about this individual's own copies are spread over one agreed set of
    /// lengths.
    pub reachable_lengths: &'a [u32],
    /// A read from another individual's DNA, or `None` for the two-term form.
    ///
    /// **On by default wherever the pre-pass emits a fraction above its floor** (spec §4.5.1,
    /// owner 2026-08-19) — the same rule the SNP/indel path follows, and for the same reason:
    /// contamination is a property of the sample, not of the marker, so a caller that corrects
    /// for it at one kind of locus and not the other is treating one number as two.
    pub contamination: Option<SsrContaminationMixture<'a>>,
}

/// **Which of the locus's reachable lengths this observation is compatible with**, as a range
/// into the ascending support.
///
/// **The three terms of the mixture must be probabilities of the same event**, and this is what
/// makes them so. A read that spanned the tract showed one length, so it allows that one; a read
/// that ran out shows a lower bound, so it allows every length at or above what it witnessed —
/// the reading spec §5.2 already gives the emission, and the one §2.1 asks for by writing the
/// junk term `λ · U(o)`, a function of the observation rather than a constant.
///
/// **Getting this wrong is not visible in a row.** Before it was written, the junk term was the
/// probability that a junk read showed *exactly* the witnessed length while the other two terms
/// were tails: measured on a locus of 31 reachable lengths, a read witnessing 8 bases got a junk
/// floor 25 times smaller than the matching tail form, so a truncated read was preferentially
/// explained as somebody else's DNA by about two orders of magnitude — and the gap widened as
/// the tract lengthened.
///
/// A complete read at a length no candidate reaches allows nothing, and the empty range says so.
fn lengths_the_observation_allows(
    observation: &SequenceObservation,
    reachable_lengths: &[u32],
) -> std::ops::Range<usize> {
    let showed = observation.bases.len() as u32;
    match observation.read_witness {
        ReadWitness::Complete => match reachable_lengths.binary_search(&showed) {
            Ok(at) => at..at + 1,
            Err(_) => 0..0,
        },
        ReadWitness::Partial { .. } => {
            reachable_lengths.partition_point(|length| *length < showed)..reachable_lengths.len()
        }
    }
}

/// **One sample's log-likelihood for every candidate genotype at one repeat tract**, written
/// into `out` in genotype-table order.
///
/// The module's own documentation carries the formula and what each term is for. This is what
/// a caller has to get right.
///
/// What the locus contributes is [`SsrLocusParameters`], which carries its own documentation for
/// each field and why none of them is derived here.
///
/// # Panics
///
/// On any mismatch between the tables handed in: a row of the wrong width; a candidate set the
/// genotype table does not cover; a context table for a different locus; a ploidy past
/// `MAX_PLOIDY_COPIES`; an outlier weight outside `(0, 1)`; an empty reachable-length support;
/// a contamination seed of a different width from that support, or one that does not sum to
/// one; or a contamination fraction outside `[0, 1)`.
///
/// **All of them are checked once here rather than per observation**, so that the accessors
/// inside the loops cannot fail in a run; a lazily checked mismatch would surface at whichever
/// locus first reached past an end, or never.
pub fn genotype_log_likelihood_row<Model: SsrEmissionModel>(
    model: &Model,
    evidence: &SsrSampleEvidence<'_>,
    locus: SsrLocusParameters<'_>,
    genotypes: &GenotypeTableView<'_>,
    out: &mut [LogProb],
    scratch: &mut SsrRowScratch<Model::Scratch>,
) {
    let SsrLocusParameters {
        candidates,
        contexts,
        outlier_weight,
        reachable_lengths,
        contamination,
    } = locus;
    let genotype_count = genotypes.genotype_count();
    let allele_count = genotypes.allele_count();
    assert_eq!(
        out.len(),
        genotype_count,
        "a row holds one entry per candidate genotype — {genotype_count}, not {}",
        out.len()
    );
    assert_eq!(
        candidates.len(),
        allele_count,
        "this locus was handed {} candidates and a genotype table over {allele_count} alleles, \
         so one of them belongs to a different locus",
        candidates.len()
    );
    assert_eq!(
        contexts.candidate_count(),
        allele_count,
        "the context table covers {} candidates and this locus is called over {allele_count}",
        contexts.candidate_count()
    );
    assert!(
        !reachable_lengths.is_empty(),
        "the outlier term is spread over the lengths the model can reach, and there is always \
         at least one — a candidate's own"
    );
    // **The support is searched, so it arrives ascending and without repeats.** The invariant
    // lives in `fill_reachable_lengths`, which sorts, and is relied on here — handing the same
    // seed a reversed slice returns 0.0 where the sorted one returns 0.335, with no panic.
    assert!(
        reachable_lengths.windows(2).all(|pair| pair[0] < pair[1]),
        "the reachable lengths are searched, so they arrive ascending and without repeats"
    );
    if let Some(contamination) = contamination {
        assert_eq!(
            contamination.contaminant_length_frequencies.len(),
            reachable_lengths.len(),
            "the contamination seed covers {} lengths and this locus reaches {} — so one of \
             them was built against a different candidate set",
            contamination.contaminant_length_frequencies.len(),
            reachable_lengths.len()
        );
        assert!(
            contamination
                .fraction_of_each_read_group
                .iter()
                .all(|view| (0.0..1.0).contains(&view.fraction)),
            "a contamination fraction is the share of a library's reads that came from \
             somebody else, so it lies in [0, 1) — not {:?}",
            contamination
                .fraction_of_each_read_group
                .iter()
                .map(|view| view.fraction)
                .collect::<Vec<_>>()
        );
        assert_eq!(
            contamination.fraction_of_each_read_group.len(),
            contexts.read_group_count(),
            "the contamination table covers {} read groups and this locus's parameters cover \
             {} — so one of them belongs to a different run",
            contamination.fraction_of_each_read_group.len(),
            contexts.read_group_count()
        );
        // **No share of a distribution is negative**, and the sum below cannot catch it: a seed
        // with one entry at −0.5 and another at +1.5 sums to one and takes the whole bracket
        // negative, which `ln` answers with `NaN` for every genotype.
        assert!(
            contamination
                .contaminant_length_frequencies
                .iter()
                .all(|share| *share >= 0.0),
            "a contamination seed is how common each length is, so no entry is negative"
        );
        // **Neither can the two shares together take everything.** Spec §4.5.1 floors
        // `1 − λ − c` positive, and the floor stays below as a guard on the arithmetic — but a
        // parameter row where the outlier weight and a contamination fraction sum past one is a
        // fit that has gone wrong, and flattening the row to a millionth silently is the
        // "confident and wrong" answer rather than the loud one.
        assert!(
            contamination
                .fraction_of_each_read_group
                .iter()
                .all(|view| outlier_weight + view.fraction < 1.0),
            "the outlier weight and a contamination fraction together take every read away from \
             this individual: {outlier_weight} plus {:?}",
            contamination
                .fraction_of_each_read_group
                .iter()
                .map(|view| view.fraction)
                .collect::<Vec<_>>()
        );
        // **The seed is a distribution and must sum to one**, and unlike its width this is
        // arithmetically load-bearing rather than a bookkeeping check: a censored read takes a
        // *suffix* of it, so a seed summing to two at a contamination fraction of 3 in 100 puts
        // 0.06 into a bracket that already holds 0.97 — a likelihood too high for every
        // genotype that observation touches, with nothing failing.
        let seed_total: f64 = contamination.contaminant_length_frequencies.iter().sum();
        assert!(
            (seed_total - 1.0).abs() < 1e-9,
            "the contamination seed is how common each length is, so it sums to one — not \
             {seed_total}"
        );
    }
    // **The outlier weight is a share of the reads, so it lives strictly inside 0 and 1**, and
    // both ends are checked because both produce a number rather than a crash. Above one,
    // `1 − λ` is negative and the bracket goes negative wherever a genotype explains the read,
    // so `ln` returns `NaN` — measured at λ = 1.5, two of three genotypes came back `NaN` and
    // the third a plausible −12.74. At exactly zero the row loses its only floor and a read
    // nothing explains takes every genotype to `−∞`, whose differences are `NaN` — which is
    // precisely what spec §4.5's junk term exists to prevent. A measurement that genuinely
    // wants no outlier term should say so rather than reach it through this argument.
    assert!(
        outlier_weight > 0.0 && outlier_weight < 1.0,
        "the outlier weight is the share of reads that came from somewhere else, so it lies \
         strictly inside 0 and 1 — not {outlier_weight}"
    );

    // **The cache is filled with `NaN` and not with zero.** An unwritten slot has to hold
    // something the row cannot mistake for a real score, and zero is exactly that mistake: a
    // slip a candidate cannot reach legitimately scores zero (spec §4.2), so zeros would make
    // *never computed* and *computed as impossible* the same value.
    scratch.prepare_emissions(evidence, allele_count, f64::NAN);
    fill_emissions(model, evidence, candidates, contexts, scratch);

    // `k / P` for every copy count a genotype can carry, shared with the SNP/indel row so the
    // two cannot disagree about what that is. It also carries the ploidy check.
    let copy_share = super::copy_shares(genotypes.ploidy());

    for slot in out.iter_mut() {
        *slot = LogProb(0.0);
    }

    // **The observations are walked in the order the caller handed them**, and this row
    // imposes none of its own — no sorting, no bucketing, no re-grouping. The caller always
    // hands it the merge's order, and that is what makes a run reproducible at any worker
    // count (spec §12 test 8).
    let counts = genotypes.genotype_allele_counts();
    for (position, observation) in evidence.observations.iter().enumerate() {
        let reads = f64::from(observation.num_obs);

        // **The junk half of the mixture, and it is a function of the observation** — spec
        // §2.1 writes it `λ · U(o)` rather than as a constant. A read no allele explains could
        // have shown any reachable length and one is as likely as another (spec §4.5), but what
        // the *observation* says is one length for a read that spanned the tract and a lower
        // bound for one that ran out, so the mass it collects differs.
        //
        // **The floor at one length is what keeps the term doing its job.** §4.5's whole reason
        // for existing is that without somewhere to put a read nothing explains, one such read
        // drives every genotype to zero — so a read whose own length is outside the candidates'
        // reach must still land somewhere, and the smallest thing a read can be is one length.
        // Here the reachable count is a **normaliser for how many lengths are in play**, not a
        // membership test; the seed below is a real distribution and does test membership,
        // which is why the two are not the same expression.
        let allowed = lengths_the_observation_allows(observation, reachable_lengths);
        let from_the_junk_distribution =
            outlier_weight * allowed.len().max(1) as f64 / reachable_lengths.len() as f64;

        // **The third term, and the reason it is a third rather than part of the second.** A
        // contaminating read shows a *plausible* length, so its distribution is peaked and its
        // term still carries the genotype through what is left over for this individual; a junk
        // read shows anything, so its term is flat and cancels between genotypes. Both are
        // properties of the observation and of the library it came from — never of the genotype
        // — so both are computed here rather than inside the genotype loop.
        let (contamination_fraction, from_the_contaminant) = match contamination {
            Some(mixture) => {
                let fraction = mixture.fraction_of(observation.read_group);
                (
                    fraction,
                    fraction * mixture.contaminant_frequency_of(observation, reachable_lengths),
                )
            }
            None => (0.0, 0.0),
        };
        // `1 − λ − c`, floored positive (spec §4.5.1).
        let from_this_individual =
            (1.0 - outlier_weight - contamination_fraction).max(MIN_SHARE_FROM_THIS_INDIVIDUAL);

        for (genotype, slot) in out.iter_mut().enumerate() {
            let carried_copies = &counts[genotype * allele_count..][..allele_count];
            // **The copy-weighted mixture over the genotype's own alleles.** A candidate no
            // copy carries is skipped rather than multiplied by zero, which is what keeps the
            // cost proportional to the ploidy rather than to the candidate count.
            let mut explained_by_this_genotype = 0.0;
            for (candidate, &copies) in carried_copies.iter().enumerate() {
                if copies == 0 {
                    continue;
                }
                explained_by_this_genotype +=
                    copy_share[copies as usize] * scratch.emission_at(position, candidate);
            }
            slot.0 += reads
                * (from_this_individual * explained_by_this_genotype
                    + from_the_junk_distribution
                    + from_the_contaminant)
                    .ln();
        }
    }
}

/// Score every `(observation, candidate)` pair into the cache, routing each observation by
/// what its reads actually witnessed.
///
/// **The witness decides the method, and nothing else does.** A read that spanned the tract
/// pins a length and goes to [`SsrEmissionModel::emission`]; a read that ran off its own end
/// proves only a lower bound and goes to [`SsrEmissionModel::censored_emission`]. Scoring a
/// partial as though it were complete mis-scores it as a *short* allele, because its bases are
/// a prefix of the truth (spec §5.1) — so the split is taken from
/// [`SsrSampleEvidence`]'s two filters, which are the single place that decides what reaches
/// the censored term, rather than re-derived from the bases here.
///
/// **The two filters together cover every observation exactly once**, which is what lets the
/// cache be filled with `NaN` and still hold no `NaN` when this returns: they enumerate the
/// same slice and split it on an exhaustive match.
fn fill_emissions<Model: SsrEmissionModel>(
    model: &Model,
    evidence: &SsrSampleEvidence<'_>,
    candidates: &[SsrCandidate<'_>],
    contexts: SsrScoringContextTable<'_>,
    scratch: &mut SsrRowScratch<Model::Scratch>,
) {
    for (position, observation) in evidence.complete_observations() {
        for (candidate, allele) in candidates.iter().enumerate() {
            let context = contexts.of(observation.read_group, candidate);
            let scored = model.emission(
                &observation.bases,
                allele,
                context,
                scratch.model_scratch_mut(),
            );
            scratch.set_emission(position, candidate, scored);
        }
    }

    for (position, observation) in evidence.partial_observations() {
        for (candidate, allele) in candidates.iter().enumerate() {
            let context = contexts.of(observation.read_group, candidate);
            let scored = model.censored_emission(
                &observation.bases,
                allele,
                context,
                scratch.model_scratch_mut(),
            );
            scratch.set_emission(position, candidate, scored);
        }
    }
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::num::NonZeroU32;

    use super::*;
    use crate::ng::alignment::StutterModel;
    use crate::ng::calling::genotype_table::GenotypeTable;
    use crate::ng::calling::likelihood::ssr_emission::fill_reachable_lengths;
    use crate::ng::calling::likelihood::ssr_emission::{
        StutterSubstitutionEmission, StutterSubstitutionScratch,
    };
    use crate::ng::calling::likelihood::stutter_rates::stutter_model_for;
    use crate::ng::locus_generation::{LocusLen, ReadWitness, SequenceObservation, SsrDetail};
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;
    use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;
    use crate::ng::types::{ErrorRate, Motif, Ploidy};

    /// The per-base substitution rate the fixtures score at — **not** a floating-point
    /// tolerance, which is what `EPSILON` would read as a few lines from `f64::EPSILON`.
    const SUBSTITUTION_RATE: f64 = 1e-3;

    fn a_motif(bases: &[u8]) -> Motif {
        Motif::new(bases).expect("a valid test motif")
    }

    fn repeats(count: u32) -> NonZeroU32 {
        NonZeroU32::new(count).expect("a test candidate always holds a repeat")
    }

    /// A slippage row at a stated level — what makes two candidates' contexts different.
    fn a_model_slipping(level: f64) -> StutterModel {
        stutter_model_for(&Slippage {
            level,
            shorter_share: 0.83,
            fall_off: 0.35,
        })
    }

    fn a_tract(motif: &[u8], repeat_count: usize) -> Vec<u8> {
        motif.repeat(repeat_count)
    }

    /// One observation: some bases, seen by `reads` reads that spanned the whole tract.
    fn spanning(bases: &[u8], reads: u32) -> SequenceObservation {
        SequenceObservation {
            bases: bases.to_vec().into_boxed_slice(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: reads,
            num_fwd: reads,
            q_sum: crate::ng::types::SummedLogError::from_nats(-10.0 * f64::from(reads)),
            mapq_sum: 60 * reads,
            mapq_sum_sq: u64::from(reads) * 3_600,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    /// The detail an STR locus carries — the motif, and flanks nothing here reads.
    fn a_detail(motif: &[u8]) -> SsrDetail {
        SsrDetail {
            motif: a_motif(motif),
            left_flank: Box::from(&b"GGGGG"[..]),
            right_flank: Box::from(&b"TTTTT"[..]),
        }
    }

    /// Everything a row needs at one locus, owned so the borrows below stay simple.
    ///
    /// **Every candidate gets its own [`StutterModel`] and its own substitution rate**, and
    /// that is not decoration. A read's chance of slipping is a property of the tract it was
    /// copied from, so candidates of different repeat counts are drawn from different strata
    /// and the row must look each one up separately (spec §4.4). A fixture that gave every
    /// candidate the same parameters could not tell that apart from a row that read candidate
    /// zero's context for all of them — and the first version of this file could not: hoisting
    /// the lookup out of the candidate loop left every row bit-identical.
    ///
    /// **The read groups differ too**, for the same reason on the other axis of the table.
    struct Fixture {
        motif: Motif,
        /// One per `(read group, candidate)`, in the table's own order.
        models: Vec<StutterModel>,
        /// One per `(read group, candidate)`, matching `models`.
        substitution_rates: Vec<f64>,
        read_groups: usize,
        candidate_bases: Vec<Vec<u8>>,
        repeat_counts: Vec<u32>,
    }

    impl Fixture {
        /// Candidates of `repeat_counts` whole copies of `motif`, in that order, scored by one
        /// read group.
        fn of(motif: &[u8], repeat_counts: &[u32]) -> Self {
            Self::of_groups(motif, repeat_counts, 1)
        }

        /// The same across `read_groups` read groups, all sharing one slippage row and one
        /// substitution rate.
        ///
        /// **For the tests that vary something else across groups.** The ordinary two-group
        /// fixture gives each group its own parameters, which is right for testing the lookup
        /// and wrong for testing anything downstream of it: the rows would differ whether or
        /// not the thing under test worked.
        fn sharing_parameters(motif: &[u8], repeat_counts: &[u32], read_groups: usize) -> Self {
            let mut fixture = Self::of_groups(motif, repeat_counts, read_groups);
            for group in 1..read_groups {
                for candidate in 0..repeat_counts.len() {
                    let at = group * repeat_counts.len() + candidate;
                    fixture.models[at] = fixture.models[candidate].clone();
                    fixture.substitution_rates[at] = fixture.substitution_rates[candidate];
                }
            }
            fixture
        }

        /// The same, across `read_groups` read groups — each with its own slippage row, as a
        /// per-chemistry fit gives them.
        fn of_groups(motif: &[u8], repeat_counts: &[u32], read_groups: usize) -> Self {
            let mut models = Vec::new();
            let mut substitution_rates = Vec::new();
            for group in 0..read_groups {
                for count in repeat_counts {
                    // Longer tracts slip more, and one lane slips more than another: two
                    // separate reasons for two contexts to differ, so a row that dropped
                    // either axis of the lookup has something to fail against.
                    //
                    // **Both are keyed by the candidate's repeat count, never by its position
                    // in the table** — that is what a stratum lookup does, and a fixture keyed
                    // by position would make permuting the candidates change the answer for a
                    // reason that has nothing to do with the row.
                    models.push(a_model_slipping(
                        0.01 + 0.004 * f64::from(*count) + 0.002 * group as f64,
                    ));
                    substitution_rates
                        .push(SUBSTITUTION_RATE * (1.0 + 0.1 * f64::from(*count) + group as f64));
                }
            }
            Self {
                motif: a_motif(motif),
                models,
                substitution_rates,
                read_groups,
                candidate_bases: repeat_counts
                    .iter()
                    .map(|count| a_tract(motif, *count as usize))
                    .collect(),
                repeat_counts: repeat_counts.to_vec(),
            }
        }

        fn candidates(&self) -> Vec<SsrCandidate<'_>> {
            self.candidate_bases
                .iter()
                .zip(&self.repeat_counts)
                .map(|(bases, count)| SsrCandidate {
                    bases,
                    repeat_count: repeats(*count),
                })
                .collect()
        }

        /// The whole `(read group, candidate)` table, in the order
        /// [`SsrScoringContextTable`] indexes it.
        fn contexts<'a>(&'a self, candidates: &[SsrCandidate<'_>]) -> Vec<SsrScoringContext<'a>> {
            let mut contexts = Vec::new();
            for group in 0..self.read_groups {
                for (candidate, allele) in candidates.iter().enumerate() {
                    let at = group * candidates.len() + candidate;
                    contexts.push(SsrScoringContext::new(
                        &self.motif,
                        &self.models[at],
                        allele,
                        ErrorRate::try_new(self.substitution_rates[at]).expect("a valid rate"),
                        [Provenance::FittedHere],
                    ));
                }
            }
            contexts
        }
    }

    /// Build the row for one set of observations, at a stated ploidy and outlier weight.
    fn score_row_at(
        fixture: &Fixture,
        observations: &[SequenceObservation],
        ploidy: u8,
        outlier_weight: f64,
    ) -> Vec<LogProb> {
        score_row_with_contamination(fixture, observations, ploidy, outlier_weight, None, &[])
    }

    /// The row at a stated outlier weight and contamination, over the locus's own reachable
    /// lengths.
    fn score_row_with_contamination(
        fixture: &Fixture,
        observations: &[SequenceObservation],
        ploidy: u8,
        outlier_weight: f64,
        seed: Option<&[f64]>,
        fractions: &[f64],
    ) -> Vec<LogProb> {
        let detail = a_detail(fixture.motif.as_bytes());
        let evidence = SsrSampleEvidence::new(observations, &detail);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let views: Vec<ContaminationView> = if fractions.is_empty() {
            vec![a_contamination_view(0.0); fixture.read_groups.max(1)]
        } else {
            fractions.iter().map(|f| a_contamination_view(*f)).collect()
        };
        let contamination = seed.map(|contaminant_length_frequencies| SsrContaminationMixture {
            fraction_of_each_read_group: &views,
            contaminant_length_frequencies,
        });
        let table = GenotypeTable::build(
            Ploidy::try_new(ploidy).expect("a valid ploidy"),
            candidates.len(),
        );
        let view = table.view();
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &StutterSubstitutionEmission,
            &evidence,
            SsrLocusParameters {
                candidates: &candidates,
                contexts: SsrScoringContextTable::new(&contexts, candidates.len()),
                outlier_weight,
                reachable_lengths: &lengths,
                contamination,
            },
            &view,
            &mut out,
            &mut scratch,
        );
        out
    }

    /// The same at the inherited outlier weight, which is what every test but two wants.
    fn score_row(
        fixture: &Fixture,
        observations: &[SequenceObservation],
        ploidy: u8,
    ) -> Vec<LogProb> {
        score_row_at(fixture, observations, ploidy, DEFAULT_OUTLIER_WEIGHT)
    }

    /// Whether no candidate at all can produce this observation — what makes a read *junk*
    /// rather than weak evidence, and the thing the cancellation test is about.
    fn every_candidate_scores_it_zero(
        fixture: &Fixture,
        observation: &SequenceObservation,
    ) -> bool {
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let mut scratch = StutterSubstitutionScratch::default();
        candidates
            .iter()
            .zip(&contexts)
            .all(|(candidate, context)| {
                StutterSubstitutionEmission.emission(
                    &observation.bases,
                    candidate,
                    context,
                    &mut scratch,
                ) == 0.0
            })
    }

    /// **Spec §12's sixth test: the junk term cancels for a read nothing explains.**
    ///
    /// A read from a paralogous tract, a chimera, a somatic length — its emission is zero under
    /// every candidate, so the bracket collapses to `λ · U`, the same number under every
    /// genotype. **What that must not do is change how far apart two genotypes are**, because
    /// that is what the caller normalises and calls with.
    ///
    /// # It cancels to one unit in the last place of the entries, not of their difference
    ///
    /// The specification asked for bitwise, and no implementation of this formula can give it:
    /// the junk term is added to each genotype's running total, and `(a + k) − (b + k)` is not
    /// `a − b` in floating point however carefully `k` is computed.
    ///
    /// **The unit matters more than the number here, and getting it wrong is how this test was
    /// first written.** Counting units in the last place *of the difference* measures the true
    /// rounding error scaled by `|entry| / |separation|`, and that ratio is set by how many junk
    /// reads there are — so the same fixture reports 16 units at 3 junk reads and 3,072 at 300,
    /// with nothing about the row having changed. Measured relative to **the entries' own
    /// magnitude**, which is where the rounding actually happens, the worst disagreement over
    /// every cell below is one `f64::EPSILON` and stays there. That is the same shape of bound
    /// `permuting_the_observations_and_the_candidates_moves_no_genotype_meaningfully` uses, and
    /// what spec §12's eighth test means by "the same relative bound as test 9".
    ///
    /// The sweep varies the junk read count for exactly that reason: a fixture with three junk
    /// reads cannot tell a robust bound from a lucky one.
    #[test]
    fn a_read_nothing_explains_moves_no_genotype_against_another() {
        let mut worst_relative = 0.0f64;
        let mut cells = 0usize;
        let mut cells_that_moved = 0usize;

        for candidate_counts in [vec![4u32, 5], vec![4, 5, 6, 7], vec![3, 6, 9]] {
            let fixture = Fixture::of(b"CA", &candidate_counts);
            let longest = *candidate_counts.iter().max().expect("a candidate");
            let explained = spanning(&a_tract(b"CA", candidate_counts[0] as usize), 7);
            let other = spanning(&a_tract(b"CA", longest as usize), 4);

            for junk_reads in [3u32, 30, 100, 300] {
                // **Eleven repeats past the *longest* candidate, not the shortest.** One past
                // the whole-repeat cutoff of a short candidate is still within the cutoff of a
                // longer one, which the distribution scores rather than refuses — so a junk
                // read built from the shortest candidate is not junk at all. The first draft of
                // this fixture made that mistake and failed for a reason that had nothing to do
                // with the property under test.
                let junk = spanning(&a_tract(b"CA", (longest + 11) as usize), junk_reads);

                // **The fixture has to be junk, or this test is about nothing.** Asserted
                // rather than reasoned about, because the reasoning is what went wrong.
                assert!(
                    every_candidate_scores_it_zero(&fixture, &junk),
                    "some candidate explains the junk read, so the term under test does not \
                     cancel"
                );

                let without = score_row(&fixture, &[explained.clone(), other.clone()], 2);
                let with = score_row(&fixture, &[explained.clone(), other.clone(), junk], 2);

                assert!(
                    with.iter().zip(&without).all(|(a, b)| a.0 < b.0),
                    "the junk read should cost every genotype something"
                );
                for first in 0..without.len() {
                    for second in 0..without.len() {
                        let moved = (with[first].0 - with[second].0)
                            - (without[first].0 - without[second].0);
                        // The entries are where the rounding happens, so they are what the
                        // error is measured against.
                        let magnitude = with[first].0.abs().max(with[second].0.abs());
                        let relative = moved.abs() / magnitude;
                        assert!(
                            relative <= f64::EPSILON,
                            "genotypes {first} and {second} moved {moved} at {junk_reads} junk \
                             reads — {relative} of an entry of magnitude {magnitude}"
                        );
                        worst_relative = worst_relative.max(relative);
                        if moved != 0.0 {
                            cells_that_moved += 1;
                        }
                        cells += 1;
                    }
                }
            }
        }

        // **The sweep has to reach cells where the two genuinely differ**, or the bound is
        // being read off arithmetic that happened to cancel and this test would pass just as
        // well against the bitwise claim it replaces.
        assert_eq!(cells, 4 * (9 + 100 + 36), "the sweep changed size");
        assert!(
            cells_that_moved > 0,
            "no cell moved at all, so this sweep cannot tell a real bound from a bitwise one"
        );
        assert!(
            worst_relative > 0.0 && worst_relative <= f64::EPSILON,
            "the worst relative disagreement was {worst_relative}"
        );
    }

    /// **Spec §12's seventh test: ploidy generality, and both of its cases.**
    ///
    /// At ploidy 4 with every read matching one allele, a genotype carrying two copies each of
    /// two alleles scores **between** the two homozygous quadruples. Where the reads are split
    /// between those alleles it scores **above both** — and that second case is the whole reason
    /// a mixed genotype is callable, so pinning only the first would leave a wrong copy
    /// weighting undetected.
    #[test]
    fn a_mixed_genotype_sits_between_the_homozygotes_and_above_them_when_the_reads_split() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let four = a_tract(b"CA", 4);
        let five = a_tract(b"CA", 5);

        let table = GenotypeTable::build(Ploidy::try_new(4).expect("a valid ploidy"), 2);
        let view = table.view();
        let counts = view.genotype_allele_counts();
        let find = |first: u32, second: u32| {
            (0..view.genotype_count())
                .find(|genotype| {
                    counts[genotype * 2] == first && counts[genotype * 2 + 1] == second
                })
                .expect("the genotype table holds every copy split")
        };
        let all_four = find(4, 0);
        let all_five = find(0, 4);
        let mixed = find(2, 2);

        // Every read matches the four-repeat allele.
        let matching = score_row(&fixture, &[spanning(&four, 10)], 4);
        assert!(
            matching[mixed].0 < matching[all_four].0 && matching[mixed].0 > matching[all_five].0,
            "with every read on one allele the mixed genotype must sit between the two \
             homozygous quadruples: {:?}",
            (
                matching[all_four].0,
                matching[mixed].0,
                matching[all_five].0
            )
        );

        // The reads split evenly between the two alleles.
        let split = score_row(&fixture, &[spanning(&four, 5), spanning(&five, 5)], 4);
        assert!(
            split[mixed].0 > split[all_four].0 && split[mixed].0 > split[all_five].0,
            "with the reads split the mixed genotype must beat both homozygous quadruples: {:?}",
            (split[all_four].0, split[mixed].0, split[all_five].0)
        );
    }

    /// **Spec §12's eighth test: the row imposes no order of its own.**
    ///
    /// It must not sort, bucket or re-group behind the caller: the caller always hands it the
    /// merge's order, and that is what makes a run reproducible at any worker count.
    ///
    /// **Not asserted bitwise, and the specification says why**: permuting the observations
    /// *is* changing the summation order, so the two rows may differ in the last bits. What is
    /// asserted is that they differ by no more than that — measured here at zero units in the
    /// last place, which is stronger than the bound and is recorded rather than relied on.
    #[test]
    fn permuting_the_observations_and_the_candidates_moves_no_genotype_meaningfully() {
        let fixture = Fixture::of(b"CA", &[4, 5, 6]);
        let observations = [
            spanning(&a_tract(b"CA", 4), 7),
            spanning(&a_tract(b"CA", 5), 3),
            spanning(&a_tract(b"CA", 6), 2),
        ];
        let forwards = score_row(&fixture, &observations, 2);

        let mut backwards_observations = observations.to_vec();
        backwards_observations.reverse();
        let backwards = score_row(&fixture, &backwards_observations, 2);

        for (genotype, (one, other)) in forwards.iter().zip(&backwards).enumerate() {
            let apart = (one.0 - other.0).abs();
            assert!(
                apart <= 2.0 * f64::EPSILON * one.0.abs().max(1.0),
                "genotype {genotype} moved {apart} nats when the observations were permuted"
            );
        }

        // **And the candidates, which the first version of this test did not do despite its
        // name.** Reversing the candidate order relabels the genotypes, so the match is made on
        // copy counts rather than by mirroring the index — at three candidates the reversal
        // maps genotype indices 0→5, 1→4, 2→2, 3→3, 4→1, 5→0, and a naive mirror compares two
        // genotypes that are twenty nats apart.
        let reversed_counts: Vec<u32> = vec![6, 5, 4];
        let reversed_fixture = Fixture::of(b"CA", &reversed_counts);
        let reversed = score_row(&reversed_fixture, &observations, 2);

        let table = GenotypeTable::build(Ploidy::try_new(2).expect("a valid ploidy"), 3);
        let view = table.view();
        let counts = view.genotype_allele_counts();
        let mut matched = 0usize;
        for forward_genotype in 0..view.genotype_count() {
            let carried = &counts[forward_genotype * 3..][..3];
            // The same genotype under the reversed candidate order carries the mirrored counts.
            let mirrored: Vec<u32> = carried.iter().rev().copied().collect();
            let reversed_genotype = (0..view.genotype_count())
                .find(|genotype| counts[genotype * 3..][..3] == mirrored[..])
                .expect("every copy split appears under either ordering");
            assert_eq!(
                forwards[forward_genotype].0.to_bits(),
                reversed[reversed_genotype].0.to_bits(),
                "genotype {carried:?} moved when the candidates were permuted"
            );
            matched += 1;
        }
        assert_eq!(matched, 6, "the genotype table changed shape");
    }

    /// **A row costs `observations × candidates` emission calls — not `× genotypes`.**
    ///
    /// Spec §8 calls that the design rather than an optimisation, and this is the only test
    /// that can tell the difference: every other one would pass just as well if the row
    /// recomputed each emission for every genotype, because the numbers would be identical.
    ///
    /// The count is instrumented rather than argued, and the fixture is chosen so the three
    /// plausible costs are three different numbers. At three observations, three candidates and
    /// a diploid — six genotypes, nine carried-allele slots — **the design costs 9 calls**,
    /// recomputing per genotype would cost 18, and recomputing per carried allele would cost
    /// 27. An earlier version of this comment gave 18 the "per carried allele" rule, which is
    /// the arithmetic for a different shape.
    #[test]
    fn one_row_scores_each_observation_against_each_candidate_exactly_once() {
        /// Counts what the row asks of it and forwards to the real model.
        struct Counting {
            inner: StutterSubstitutionEmission,
            complete: Cell<usize>,
            censored: Cell<usize>,
        }

        impl SsrEmissionModel for Counting {
            type Scratch = StutterSubstitutionScratch;

            fn emission(
                &self,
                observation: &[u8],
                candidate: &SsrCandidate<'_>,
                context: &SsrScoringContext<'_>,
                scratch: &mut Self::Scratch,
            ) -> f64 {
                self.complete.set(self.complete.get() + 1);
                self.inner
                    .emission(observation, candidate, context, scratch)
            }

            fn censored_emission(
                &self,
                witnessed_prefix: &[u8],
                candidate: &SsrCandidate<'_>,
                context: &SsrScoringContext<'_>,
                scratch: &mut Self::Scratch,
            ) -> f64 {
                self.censored.set(self.censored.get() + 1);
                self.inner
                    .censored_emission(witnessed_prefix, candidate, context, scratch)
            }
        }

        let fixture = Fixture::of(b"CA", &[4, 5, 6]);
        let observations = [
            spanning(&a_tract(b"CA", 4), 7),
            spanning(&a_tract(b"CA", 5), 3),
            spanning(&a_tract(b"CA", 6), 2),
        ];
        let detail = a_detail(b"CA");
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let table = GenotypeTable::build(Ploidy::try_new(2).expect("a valid ploidy"), 3);
        let view = table.view();
        assert_eq!(
            view.genotype_count(),
            6,
            "the fixture must hold more genotypes than candidates, or it cannot tell the two \
             costs apart"
        );

        let model = Counting {
            inner: StutterSubstitutionEmission,
            complete: Cell::new(0),
            censored: Cell::new(0),
        };
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &model,
            &evidence,
            SsrLocusParameters {
                candidates: &candidates,
                contexts: SsrScoringContextTable::new(&contexts, candidates.len()),
                outlier_weight: DEFAULT_OUTLIER_WEIGHT,
                reachable_lengths: &lengths,
                contamination: None,
            },
            &view,
            &mut out,
            &mut scratch,
        );

        assert_eq!(
            model.complete.get(),
            observations.len() * candidates.len(),
            "the row scored each (observation, candidate) pair {} times",
            model.complete.get()
        );
        assert_eq!(model.censored.get(), 0, "no observation here is censored");
    }

    /// **Spec §12's seventh test, first half: a biallelic diploid reproduced by hand.**
    ///
    /// The other ploidy test pins *orderings* — which genotype beats which — and every one of
    /// those survives a copy weighting that is wrong by a constant factor. This one recomputes
    /// the whole formula outside the row, from the two emissions and the seven numbers spec
    /// §2.1 names, and requires the answer to the bit.
    ///
    /// **It is the test that pins `k / P` itself.** Without it, deleting the division by the
    /// ploidy — so the weights sum to the ploidy instead of to one — passes every other test in
    /// this file while moving every entry of the row by about 7 nats, a factor of a thousand in
    /// likelihood.
    #[test]
    fn a_biallelic_diploid_row_matches_the_formula_computed_by_hand() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);

        // **One of the three is a read that ran out**, which is what makes this the oracle for
        // the junk term's shape as well as for the copy weights: its `U(o)` is a tail where a
        // whole read's is one length, and nothing else in this file can see the difference.
        let mut truncated = spanning(&a_tract(b"CA", 3), 4);
        truncated.read_witness =
            ReadWitness::from_left(6, LocusLen::from_positions(10)).expect("a partial witness");
        let observations = [
            spanning(&a_tract(b"CA", 4), 6),
            spanning(&a_tract(b"CA", 5), 2),
            truncated,
        ];
        // The emissions, taken straight from the model — the row's only input this test does
        // not recompute, because it is the seam's business and not the row's.
        let mut model_scratch = StutterSubstitutionScratch::default();
        let emission = |observation: usize, candidate: usize, scratch: &mut _| -> f64 {
            let seen = &observations[observation];
            match seen.read_witness {
                ReadWitness::Complete => StutterSubstitutionEmission.emission(
                    &seen.bases,
                    &candidates[candidate],
                    &contexts[candidate],
                    scratch,
                ),
                ReadWitness::Partial { .. } => StutterSubstitutionEmission.censored_emission(
                    &seen.bases,
                    &candidates[candidate],
                    &contexts[candidate],
                    scratch,
                ),
            }
        };
        let scored: Vec<Vec<f64>> = (0..observations.len())
            .map(|observation| {
                (0..candidates.len())
                    .map(|candidate| emission(observation, candidate, &mut model_scratch))
                    .collect()
            })
            .collect();

        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let own = 1.0 - DEFAULT_OUTLIER_WEIGHT;
        // `λ · U(o)` computed here, from the specification's own definition of `U(o)` — the
        // uniform distribution's mass over the lengths this observation allows. A whole read
        // allows one; a read that ran out allows every length at or above what it witnessed.
        let junk_of = |observation: &SequenceObservation| -> f64 {
            let showed = observation.bases.len() as u32;
            let allowed = match observation.read_witness {
                ReadWitness::Complete => usize::from(lengths.contains(&showed)),
                ReadWitness::Partial { .. } => {
                    lengths.iter().filter(|length| **length >= showed).count()
                }
            };
            DEFAULT_OUTLIER_WEIGHT * allowed.max(1) as f64 / lengths.len() as f64
        };
        // The truncated read must actually allow more than one length, or this fixture cannot
        // tell a tail from a point mass.
        assert!(
            junk_of(&observations[2]) > junk_of(&observations[0]) * 1.5,
            "the truncated read's junk mass is {} against a whole read's {}",
            junk_of(&observations[2]),
            junk_of(&observations[0])
        );
        // The three diploid genotypes over two alleles, in genotype-table order: (2,0), (1,1),
        // (0,2). `k / P` is 1.0, 0.5 and 0.0 for the first allele, and the complement for the
        // second.
        let by_hand: Vec<f64> = [(1.0, 0.0), (0.5, 0.5), (0.0, 1.0)]
            .iter()
            .map(|(first, second)| {
                observations
                    .iter()
                    .enumerate()
                    .map(|(observation, entry)| {
                        let explained =
                            first * scored[observation][0] + second * scored[observation][1];
                        f64::from(entry.num_obs) * (own * explained + junk_of(entry)).ln()
                    })
                    .sum()
            })
            .collect();

        let row = score_row(&fixture, &observations, 2);
        assert_eq!(row.len(), by_hand.len(), "the genotype table changed shape");
        for (genotype, (slot, expected)) in row.iter().zip(&by_hand).enumerate() {
            assert_eq!(
                slot.0.to_bits(),
                expected.to_bits(),
                "genotype {genotype}: the row gave {} and the formula {expected}",
                slot.0
            );
        }

        // And the copy weights really are what the fixture claims: a row whose weights summed
        // to the ploidy rather than to one would be about seven nats away from this.
        assert!(
            (row[0].0 - (-11.465_325)).abs() < 1e-5,
            "the fixture moved: {}",
            row[0].0
        );
    }

    /// **A candidate is scored against its own stratum's parameters, not candidate zero's.**
    ///
    /// Spec §4.4 is explicit that the context is built per `(read group, candidate)` and that
    /// the lookup may not be hoisted out of the candidate loop: candidates of different repeat
    /// counts slip at measurably different rates. **Hoisting it is a one-character edit that
    /// nothing else in this file can see** — every other fixture would give the same answer,
    /// because they would all be reading identical parameters.
    ///
    /// The same on the other axis: two read groups with different slippage rows must score the
    /// same bases differently.
    #[test]
    fn each_candidate_and_read_group_is_scored_against_its_own_parameters() {
        let fixture = Fixture::of_groups(b"CA", &[4, 6], 2);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let table = SsrScoringContextTable::new(&contexts, candidates.len());
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let detail = a_detail(b"CA");

        // The table's four entries are four different parameter sets, which is what makes the
        // rows below able to disagree at all.
        assert_eq!(table.read_group_count(), 2);
        for group in 0..2u32 {
            let first = table.of(ReadGroupId(group), 0);
            let second = table.of(ReadGroupId(group), 1);
            assert!(
                first.stutter.same_length_share() != second.stutter.same_length_share(),
                "the two candidates of read group {group} share a slippage row, so this test \
                 cannot see a hoisted lookup"
            );
        }
        assert!(
            table.of(ReadGroupId(0), 0).stutter.same_length_share()
                != table.of(ReadGroupId(1), 0).stutter.same_length_share(),
            "the two read groups share a slippage row"
        );

        // A read of the same bases from the two read groups must be scored differently.
        let bases = a_tract(b"CA", 5);
        let from_first = {
            let mut observation = spanning(&bases, 5);
            observation.read_group = ReadGroupId(0);
            observation
        };
        let from_second = {
            let mut observation = spanning(&bases, 5);
            observation.read_group = ReadGroupId(1);
            observation
        };

        let score_with = |observation: &SequenceObservation| {
            let held = [observation.clone()];
            let evidence = SsrSampleEvidence::new(&held, &detail);
            let genotypes = GenotypeTable::build(
                Ploidy::try_new(2).expect("a valid ploidy"),
                candidates.len(),
            );
            let view = genotypes.view();
            let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
            let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
            genotype_log_likelihood_row(
                &StutterSubstitutionEmission,
                &evidence,
                SsrLocusParameters {
                    candidates: &candidates,
                    contexts: table,
                    outlier_weight: DEFAULT_OUTLIER_WEIGHT,
                    reachable_lengths: &lengths,
                    contamination: None,
                },
                &view,
                &mut out,
                &mut scratch,
            );
            out
        };

        let first = score_with(&from_first);
        let second = score_with(&from_second);
        assert!(
            first.iter().zip(&second).any(|(a, b)| a.0 != b.0),
            "the two read groups scored the same read identically, so the row is not reading \
             the read group's own parameters"
        );
    }

    /// **A read that ran out inside the tract reaches the censored term, and is not scored as a
    /// short allele.**
    ///
    /// Without a partial observation in a fixture, three separate defects pass every other test
    /// in this file: routing the partial loop to `emission`, indexing the emission cache by a
    /// dense counter inside the filtered loop, and any change to how a witness is read.
    ///
    /// **The bases are deliberately a prefix of the longer candidate**, which is the case the
    /// whole censored term exists for: scored as a complete read they say *this sample carries
    /// the short allele*, and scored as a lower bound they say *the tract is at least this
    /// long*, which the longer candidate satisfies outright.
    #[test]
    fn a_read_that_ran_out_is_scored_as_a_lower_bound_and_not_as_a_short_allele() {
        let fixture = Fixture::of(b"CA", &[4, 8]);
        let detail = a_detail(b"CA");
        let witnessed = a_tract(b"CA", 4);

        let mut partial = spanning(&witnessed, 9);
        partial.read_witness = ReadWitness::from_left(8, LocusLen::from_positions(16))
            .expect("a partial witness of eight of sixteen positions");
        assert!(
            matches!(partial.read_witness, ReadWitness::Partial { .. }),
            "the fixture must actually be partial"
        );
        let complete = spanning(&witnessed, 9);

        // The cache is indexed by position in the *whole* observation slice, so the partial is
        // put first and a complete observation above it: a dense counter inside the filtered
        // loop addresses the wrong row for that one, and nothing else here would notice.
        let observations = [partial, spanning(&a_tract(b"CA", 8), 4)];
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        assert_eq!(evidence.partial_observations().count(), 1);
        assert_eq!(evidence.complete_observations().count(), 1);

        let as_lower_bound = score_row(&fixture, &observations, 2);
        let as_short_allele = score_row(&fixture, &[complete, spanning(&a_tract(b"CA", 8), 4)], 2);

        // Genotype order over two alleles at ploidy 2: (2,0), (1,1), (0,2).
        let eight_eight = 2;
        let four_four = 0;
        assert!(
            as_lower_bound[eight_eight].0 > as_short_allele[eight_eight].0,
            "read as a lower bound, eight bases of tract should cost the eight-repeat \
             homozygote less than reading it as a whole four-repeat allele: {} against {}",
            as_lower_bound[eight_eight].0,
            as_short_allele[eight_eight].0
        );
        // **What separates the two readings is which allele the read is evidence *for*, not
        // the absolute score.** A lower bound is never less likely than the complete read of the
        // same bases under *any* candidate — that is the censored term's own invariant — so both
        // entries rise. The thing that must move is how far apart the two homozygotes are.
        let separated_as_a_short_allele =
            as_short_allele[four_four].0 - as_short_allele[eight_eight].0;
        let separated_as_a_lower_bound =
            as_lower_bound[four_four].0 - as_lower_bound[eight_eight].0;
        assert!(
            separated_as_a_short_allele > separated_as_a_lower_bound,
            "read as a whole allele the bases should favour the four-repeat homozygote more \
             than the same bases read as a lower bound do: {separated_as_a_short_allele} \
             against {separated_as_a_lower_bound}"
        );
    }

    /// **The row refuses an outlier weight that is not a share of the reads.**
    ///
    /// Both ends produce a number rather than a crash, which is why both are checked: above one
    /// the explained half of the mixture goes negative and `ln` returns `NaN`, and at exactly
    /// zero a read nothing explains takes every genotype to `−∞`, whose differences are `NaN` —
    /// the collapse spec §4.5's junk term exists to prevent.
    #[test]
    #[should_panic(expected = "the outlier weight is the share of reads")]
    fn an_outlier_weight_above_one_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row_at(&fixture, &observations, 2, 1.5);
    }

    #[test]
    #[should_panic(expected = "the outlier weight is the share of reads")]
    fn an_outlier_weight_of_zero_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row_at(&fixture, &observations, 2, 0.0);
    }

    /// **A locus reaching no lengths at all is refused**, because the outlier weight is spread
    /// over them and there is always at least one — a candidate's own.
    ///
    /// Built by hand rather than through the helper, because
    /// [`fill_reachable_lengths`] cannot produce an empty support from a real candidate set:
    /// this is the caller who passed a buffer they forgot to fill.
    #[test]
    #[should_panic(expected = "the outlier term is spread over the lengths")]
    fn a_locus_reaching_no_lengths_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        let detail = a_detail(b"CA");
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let table = GenotypeTable::build(Ploidy::try_new(2).expect("a valid ploidy"), 2);
        let view = table.view();
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &StutterSubstitutionEmission,
            &evidence,
            SsrLocusParameters {
                candidates: &candidates,
                contexts: SsrScoringContextTable::new(&contexts, candidates.len()),
                outlier_weight: DEFAULT_OUTLIER_WEIGHT,
                reachable_lengths: &[],
                contamination: None,
            },
            &view,
            &mut out,
            &mut scratch,
        );
    }

    /// **A candidate set that does not match the genotype table is refused**, rather than
    /// scoring a row over one locus's alleles and another's genotypes.
    #[test]
    #[should_panic(expected = "belongs to a different locus")]
    fn a_candidate_set_the_genotype_table_does_not_cover_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5, 6]);
        let detail = a_detail(b"CA");
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        let evidence = SsrSampleEvidence::new(&observations, &detail);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        // A genotype table over two alleles, against three candidates.
        let table = GenotypeTable::build(Ploidy::try_new(2).expect("a valid ploidy"), 2);
        let view = table.view();
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut scratch = SsrRowScratch::<StutterSubstitutionScratch>::default();
        genotype_log_likelihood_row(
            &StutterSubstitutionEmission,
            &evidence,
            SsrLocusParameters {
                candidates: &candidates,
                contexts: SsrScoringContextTable::new(&contexts, candidates.len()),
                outlier_weight: DEFAULT_OUTLIER_WEIGHT,
                reachable_lengths: &lengths,
                contamination: None,
            },
            &view,
            &mut out,
            &mut scratch,
        );
    }

    /// **A sample with no observations at this tract leaves every genotype at zero**, without a
    /// branch: the sum is empty, and the prior is what decides. Pinned because an empty row is
    /// the shape a caller is most likely to hand in without meaning to.
    #[test]
    fn a_sample_with_no_reads_leaves_every_genotype_at_zero() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let row = score_row(&fixture, &[], 2);
        assert!(
            row.iter().all(|slot| slot.0 == 0.0),
            "an empty row should be all zeros, not {row:?}"
        );
    }

    /// **The lengths a locus reaches come from the candidates and the cutoffs, and from nothing
    /// else** — no observation, no sample count, and none of the fitted rates.
    ///
    /// This is the repair of production's cohort-wide `D` (spec §4.5), so the property that
    /// matters is a negative one: the answer does not move when the reads move. A count that
    /// asked what samples showed would give one number for a single sample and a smaller one
    /// for a panel, which is how production's junk floor came to be ten times lower at 63
    /// accessions than at one.
    #[test]
    fn the_reachable_lengths_depend_on_the_candidates_and_the_cutoffs_alone() {
        let fixture = Fixture::of(b"CA", &[4, 6]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);

        // Ascending, no repeats, and the two candidates' own lengths are among them.
        assert!(lengths.windows(2).all(|pair| pair[0] < pair[1]));
        assert!(lengths.contains(&8) && lengths.contains(&12));

        // **A reused buffer is cleared**, or a second locus inherits the first one's lengths.
        // Checked against a fresh buffer rather than against a length that happens to be absent
        // — the two loci below share several lengths, so "does it still hold 12" proves nothing.
        let other = Fixture::of(b"CAG", &[3]);
        let other_candidates = other.candidates();
        fill_reachable_lengths(&other_candidates, &other.motif, &mut lengths);
        let mut into_a_fresh_buffer = Vec::new();
        fill_reachable_lengths(&other_candidates, &other.motif, &mut into_a_fresh_buffer);
        assert_eq!(
            lengths, into_a_fresh_buffer,
            "the reused buffer kept something from the previous locus"
        );

        // And nothing from the first locus that the second cannot reach survived.
        let first_only: Vec<u32> = {
            let mut first = Vec::new();
            fill_reachable_lengths(&candidates, &fixture.motif, &mut first);
            first
                .into_iter()
                .filter(|length| !into_a_fresh_buffer.contains(length))
                .collect()
        };
        assert!(
            !first_only.is_empty(),
            "the two fixtures reach the same lengths, so this check cannot see a stale buffer"
        );
        assert!(first_only.iter().all(|length| !lengths.contains(length)));
    }

    /// **The junk floor is a property of the locus, not of what the cohort showed** — spec
    /// §4.5's repair of production's cohort-wide `D`.
    ///
    /// **Pinned by its size rather than by re-running the row**, because the obvious test does
    /// not work: calling a pure function twice with the same arguments cannot fail, and the
    /// version of this test that did so killed none of seventeen mutations. What has teeth is
    /// the number itself. This sample shows two distinct sequences; production would spread the
    /// outlier weight over those two and give a floor of 0.005, and a 63-accession panel showing
    /// twenty would get 0.0005. Here the locus reaches 31 lengths whatever anybody showed, so
    /// the floor is 0.01/31 either way.
    #[test]
    fn the_junk_floor_comes_from_the_locus_and_not_from_what_was_seen() {
        let fixture = Fixture::of(b"CA", &[4, 6]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        assert_eq!(
            lengths.len(),
            31,
            "the fixture's support changed size, and the floor below is read off it"
        );

        // Two distinct sequences, which is what production would have counted.
        let observations = [
            spanning(&a_tract(b"CA", 4), 9),
            spanning(&a_tract(b"CA", 6), 7),
        ];
        let row = score_row(&fixture, &observations, 2);

        // A read nothing explains, scored against the same locus: its whole bracket is the junk
        // floor, so the row reads it straight back.
        let junk = spanning(&a_tract(b"CA", 6 + 11), 1);
        assert!(every_candidate_scores_it_zero(&fixture, &junk));
        let with_junk = score_row(
            &fixture,
            &[observations[0].clone(), observations[1].clone(), junk],
            2,
        );

        let floor = (with_junk[0].0 - row[0].0).exp();
        let from_the_locus = DEFAULT_OUTLIER_WEIGHT / 31.0;
        assert!(
            (floor - from_the_locus).abs() < 1e-12,
            "the junk floor is {floor}, not the locus's own {from_the_locus}"
        );
        // What production would have given the same sample, and the same panel.
        assert!(
            (floor - DEFAULT_OUTLIER_WEIGHT / 2.0).abs() > 1e-3,
            "the floor is production's two-sequence answer, so nothing was repaired"
        );
    }

    /// **At a contamination fraction of zero the three-term form is the two-term form**, which
    /// is what lets the term default on: it pins that a clean sample is untouched.
    ///
    /// Bit for bit here, and that is arithmetic rather than luck: `1 − λ − 0` is `1 − λ`
    /// exactly, and adding `0 · seed` adds a zero. It is the SNP/indel path's identity that
    /// needs a tolerance, because there an `exp`/`log` round trip sits between the two forms.
    #[test]
    fn a_contamination_fraction_of_zero_is_the_two_term_form() {
        let fixture = Fixture::of(b"CA", &[4, 5, 6]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let seed = a_seed_peaked_at(&lengths, 10);

        let observations = [
            spanning(&a_tract(b"CA", 4), 6),
            spanning(&a_tract(b"CA", 5), 3),
        ];
        let two_term = score_row(&fixture, &observations, 2);
        let three_term = score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[0.0],
        );

        for (genotype, (without, with)) in two_term.iter().zip(&three_term).enumerate() {
            assert_eq!(
                without.0.to_bits(),
                with.0.to_bits(),
                "genotype {genotype} moved at a contamination fraction of zero"
            );
        }
    }

    /// **What the term is for: a handful of reads carrying an allele this sample does not have
    /// stops being forced to choose between error and heterozygosity.**
    ///
    /// Spec §4.5.1: today a contaminating read has to be explained as slippage where the model
    /// can reach it — inflating the apparent slip rate — or falls to the outlier floor where it
    /// cannot, and the first is how a contaminated sample gets called heterozygous for its
    /// contaminant's allele. With the term on, those reads have a third explanation carrying a
    /// measured weight, so the heterozygote's advantage over the true homozygote shrinks.
    #[test]
    fn contamination_softens_the_pull_of_a_few_reads_carrying_another_allele() {
        let fixture = Fixture::of(b"CA", &[4, 6]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        // The contaminant's DNA is common stuff: the seed is peaked at the six-repeat length,
        // which is what the intruding reads show.
        let seed = a_seed_peaked_at(&lengths, 12);

        // Mostly this sample's own four-repeat allele, with a few reads of the other.
        let observations = [
            spanning(&a_tract(b"CA", 4), 30),
            spanning(&a_tract(b"CA", 6), 3),
        ];

        // Genotype order over two alleles at ploidy 2: (2,0), (1,1), (0,2).
        let homozygous = 0;
        let heterozygous = 1;

        let clean = score_row(&fixture, &observations, 2);
        let contaminated = score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[0.03],
        );

        let pull_when_clean = clean[heterozygous].0 - clean[homozygous].0;
        let pull_when_contaminated = contaminated[heterozygous].0 - contaminated[homozygous].0;
        assert!(
            pull_when_contaminated < pull_when_clean,
            "the heterozygote's advantage should shrink once the intruding reads have a third \
             explanation: {pull_when_contaminated} against {pull_when_clean}"
        );
    }

    /// **A read that ran out inside the tract gets the seed's mass at or above what it
    /// witnessed**, not the mass at that exact length — the same lower-bound reading spec §5.2
    /// gives the emission, and the resolution of the open question the architecture records.
    ///
    /// Without it a truncated read would be charged as a contaminant of the *short* allele,
    /// which is the same trap §5.1 names for the emission.
    #[test]
    fn a_censored_read_takes_the_seeds_mass_at_or_above_what_it_witnessed() {
        let fixture = Fixture::of(b"CA", &[4, 8]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);

        // All the seed's mass sits well above what the read witnessed, so a point lookup would
        // give the read nothing and the survival reading gives it everything.
        let mut seed = vec![0.0; lengths.len()];
        let far_above = lengths
            .iter()
            .position(|length| *length >= 16)
            .expect("a length");
        seed[far_above] = 1.0;

        let witnessed = a_tract(b"CA", 4);
        let mut partial = spanning(&witnessed, 5);
        partial.read_witness =
            ReadWitness::from_left(8, LocusLen::from_positions(16)).expect("a partial witness");

        let fractions = [a_contamination_view(0.03)];
        let contamination = SsrContaminationMixture {
            fraction_of_each_read_group: &fractions,
            contaminant_length_frequencies: &seed,
        };
        assert_eq!(
            contamination.contaminant_frequency_of(&partial, &lengths),
            1.0,
            "a lower bound should collect every length at or above it"
        );

        // The same bases seen whole collect only their own length's share, which here is zero.
        let complete = spanning(&witnessed, 5);
        assert_eq!(
            contamination.contaminant_frequency_of(&complete, &lengths),
            0.0,
            "a complete read should collect its own length's share and nothing else"
        );
    }

    /// **A contamination seed built against another locus's candidates is refused**, rather than
    /// spreading one locus's lengths over another's support.
    #[test]
    #[should_panic(expected = "built against a different candidate set")]
    fn a_seed_of_the_wrong_width_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&[0.5, 0.5]),
            &[0.02],
        );
    }

    /// **A contamination fraction outside `[0, 1)` is refused.** It is a share of a library's
    /// reads, and at one or above there is nothing left for the sample itself.
    #[test]
    #[should_panic(expected = "so it lies in [0, 1)")]
    fn a_contamination_fraction_of_one_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let seed = a_seed_peaked_at(&lengths, 8);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[1.0],
        );
    }

    /// A contamination estimate at a stated fraction — the shape the pre-pass hands over.
    fn a_contamination_view(fraction: f64) -> ContaminationView {
        ContaminationView {
            fraction,
            markers_with_reads: 1_000,
            reads_on_markers: 30_000,
            source: ContaminationSource::ThisReadGroupsReads,
        }
    }

    /// A seed distribution over `lengths` with most of its mass on one of them — the shape a
    /// stratum's fitted length spectrum produces where most of its chromosomes carry one length
    /// (`doc/devel/ng/spec/population_diversity.md` §4.2, which replaced the geometric decay
    /// away from the cohort's modal length this fixture was written against).
    fn a_seed_peaked_at(lengths: &[u32], peak: u32) -> Vec<f64> {
        let at = lengths
            .iter()
            .position(|length| *length == peak)
            .expect("the peak must be a length this locus reaches");
        let mut seed = vec![0.0; lengths.len()];
        for (index, share) in seed.iter_mut().enumerate() {
            let away = index.abs_diff(at) as i32;
            *share = 0.5f64.powi(away);
        }
        let total: f64 = seed.iter().sum();
        for share in &mut seed {
            *share /= total;
        }
        seed
    }

    /// **A seed that is not a distribution is refused**, and this one is arithmetically
    /// load-bearing rather than bookkeeping.
    ///
    /// A read that ran out takes a **suffix** of the seed, so a seed summing to two at a
    /// contamination fraction of 3 in 100 puts 0.06 into a bracket that already holds 0.97 —
    /// measured, that moved the homozygous-reference genotype from −0.4740 to −0.1643 nats,
    /// about 1.3 on the Phred scale, with nothing failing. The producer that would guarantee
    /// the sum does not exist yet, which is exactly when a check earns its place.
    #[test]
    #[should_panic(expected = "so it sums to one")]
    fn a_seed_that_is_not_a_distribution_is_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let mut seed = a_seed_peaked_at(&lengths, 8);
        // Twice the mass it should carry.
        for share in &mut seed {
            *share *= 2.0;
        }
        let observations = [spanning(&a_tract(b"CA", 4), 5)];
        score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[0.03],
        );
    }

    /// **The outlier weight and a contamination fraction that together take every read are
    /// refused**, rather than silently flattening the row.
    ///
    /// The two come from different fits and neither knows about the other, so nothing but this
    /// stops them summing past one. Spec §4.5.1 asks for `1 − λ − c` to be floored positive and
    /// it still is — but the floor is a guard on the arithmetic, not a policy: a row flattened
    /// to a millionth by a parameter set that has gone wrong is the confident-and-wrong answer,
    /// and this is the loud one.
    #[test]
    #[should_panic(expected = "together take every read away from this individual")]
    fn the_two_shares_together_taking_every_read_are_refused() {
        let fixture = Fixture::of(b"CA", &[4, 5]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let seed = a_seed_peaked_at(&lengths, 8);
        let observations = [spanning(&a_tract(b"CA", 4), 5)];

        // A contamination fraction that leaves the outlier weight nothing to sit inside.
        score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[0.995],
        );
    }

    /// **The seed's two readings, pinned by value.** Neither branch had a test that could see an
    /// off-by-one: the censored fixture put zero mass at the witnessed length, so moving the
    /// boundary changed nothing it asserted.
    ///
    /// Measured against a seed whose shares are known: a complete read collects **its own
    /// length's share and nothing else**, and a censored one collects **every share at or above
    /// what it witnessed, its own included**. Shifting either boundary by one halves or doubles
    /// the answer, which is what the numbers below say.
    #[test]
    fn the_seed_is_read_at_one_length_for_a_whole_read_and_from_it_upward_for_a_truncated_one() {
        let fixture = Fixture::of(b"CA", &[4, 6]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);

        // A seed with a different, known share at every length, so any boundary shift moves the
        // answer: halving powers of two, normalised.
        let at_eight = lengths
            .iter()
            .position(|length| *length == 8)
            .expect("eight bases is reachable");
        let mut seed = vec![0.0; lengths.len()];
        for (index, share) in seed.iter_mut().enumerate() {
            *share = 0.5f64.powi(index as i32);
        }
        let total: f64 = seed.iter().sum();
        for share in &mut seed {
            *share /= total;
        }
        let fractions = [a_contamination_view(0.03)];
        let contamination = SsrContaminationMixture {
            fraction_of_each_read_group: &fractions,
            contaminant_length_frequencies: &seed,
        };

        let witnessed = a_tract(b"CA", 4);
        let complete = spanning(&witnessed, 5);
        assert_eq!(complete.bases.len(), 8);
        let mut truncated = spanning(&witnessed, 5);
        truncated.read_witness =
            ReadWitness::from_left(8, LocusLen::from_positions(16)).expect("a partial witness");

        // The whole read collects one entry.
        let one_length = contamination.contaminant_frequency_of(&complete, &lengths);
        assert!(
            (one_length - seed[at_eight]).abs() < 1e-15,
            "a whole read collected {one_length}, not its own length's {}",
            seed[at_eight]
        );

        // The truncated read collects that entry and every one above it — which for a halving
        // seed is very nearly twice as much, and would be exactly half as much if the boundary
        // excluded the witnessed length itself.
        let from_there_up = contamination.contaminant_frequency_of(&truncated, &lengths);
        let expected: f64 = seed[at_eight..].iter().sum();
        assert!(
            (from_there_up - expected).abs() < 1e-15,
            "a truncated read collected {from_there_up}, not {expected}"
        );
        assert!(
            from_there_up > one_length * 1.9,
            "the two readings are too close for this fixture to see a boundary shift: \
             {from_there_up} against {one_length}"
        );
        // And the boundary is inclusive: excluding the witnessed length would drop the largest
        // share in the tail.
        let excluding_its_own = expected - seed[at_eight];
        assert!(
            (from_there_up - excluding_its_own).abs() > 0.4 * from_there_up,
            "excluding the witnessed length would barely move the answer, so this cannot see it"
        );
    }

    /// **A contaminating read is not a junk read, and the row must not be able to pretend it
    /// is.** Spec §4.5.1 names the wrong answer outright: adding `c` to `λ` and keeping one flat
    /// distribution gives a contaminant's reads the junk treatment, which is what the model does
    /// today and what the term exists to stop.
    ///
    /// **Nothing else in this file can tell the two apart.** Folding `c` into `λ` leaves every
    /// other test green while moving a genotype by nats, and so does dropping the seed. What
    /// separates them is that the seed is *peaked*: a read at the length the contaminant carries
    /// gets far more from `c · seed(o)` than the same weight spread flat would give it.
    #[test]
    fn the_contamination_term_is_peaked_and_the_junk_term_is_flat() {
        let fixture = Fixture::of(b"CA", &[4, 6]);
        let candidates = fixture.candidates();
        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let seed = a_seed_peaked_at(&lengths, 12);
        let peak = lengths
            .iter()
            .position(|length| *length == 12)
            .expect("reachable");
        // Peaked means *against the flat alternative*, which is the comparison this test makes:
        // a tenth of the mass on one of 31 lengths is ten times what spreading it flat gives.
        assert!(
            seed[peak] > 10.0 / lengths.len() as f64,
            "the seed must be peaked, or this test cannot tell it from a flat one: {} against a \
             flat {}",
            seed[peak],
            1.0 / lengths.len() as f64
        );

        let observations = [
            spanning(&a_tract(b"CA", 4), 30),
            spanning(&a_tract(b"CA", 6), 3),
        ];
        let fraction = 0.03;

        let three_term = score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[fraction],
        );
        // The wrong answer spec §4.5.1 names: one flat distribution carrying both weights.
        let folded_together = score_row_at(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT + fraction,
        );

        assert!(
            three_term
                .iter()
                .zip(&folded_together)
                .any(|(peaked, flat)| (peaked.0 - flat.0).abs() > 1.0),
            "the three-term mixture and `λ + c` under one flat distribution agree, so the term \
             is doing nothing a wider outlier weight would not: {three_term:?} against \
             {folded_together:?}"
        );

        // And it is the *seed* doing it: with the same fraction and a flat seed the two agree
        // much more closely.
        let flat_seed = vec![1.0 / lengths.len() as f64; lengths.len()];
        let flat_three_term = score_row_with_contamination(
            &fixture,
            &observations,
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&flat_seed),
            &[fraction],
        );
        for (genotype, (flat, folded)) in flat_three_term.iter().zip(&folded_together).enumerate() {
            assert!(
                (flat.0 - folded.0).abs() < 1e-9,
                "genotype {genotype}: a flat seed should reproduce `λ + c` exactly, {} against {}",
                flat.0,
                folded.0
            );
        }
    }

    /// **Each read group is charged its own contamination fraction.** A library sequenced clean
    /// and one carrying another individual's DNA must not be scored alike, and nothing else here
    /// varies the fraction across groups.
    ///
    /// **The two groups are given identical slippage and substitution parameters**, which is the
    /// trap: the ordinary two-group fixture gives them different ones, so the rows differ for a
    /// reason that has nothing to do with contamination and the test passes under the mutant.
    #[test]
    fn each_read_group_is_charged_its_own_contamination_fraction() {
        let fixture = Fixture::sharing_parameters(b"CA", &[4, 6], 2);
        let candidates = fixture.candidates();
        let contexts = fixture.contexts(&candidates);
        let table = SsrScoringContextTable::new(&contexts, candidates.len());
        assert_eq!(
            table.of(ReadGroupId(0), 0).stutter.same_length_share(),
            table.of(ReadGroupId(1), 0).stutter.same_length_share(),
            "the two groups must share their parameters, or this test cannot see the fraction"
        );

        let mut lengths = Vec::new();
        fill_reachable_lengths(&candidates, &fixture.motif, &mut lengths);
        let seed = a_seed_peaked_at(&lengths, 12);

        let mut from_the_clean_group = spanning(&a_tract(b"CA", 6), 6);
        from_the_clean_group.read_group = ReadGroupId(0);
        let mut from_the_dirty_group = spanning(&a_tract(b"CA", 6), 6);
        from_the_dirty_group.read_group = ReadGroupId(1);

        let clean = score_row_with_contamination(
            &fixture,
            &[spanning(&a_tract(b"CA", 4), 20), from_the_clean_group],
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[0.0, 0.20],
        );
        let dirty = score_row_with_contamination(
            &fixture,
            &[spanning(&a_tract(b"CA", 4), 20), from_the_dirty_group],
            2,
            DEFAULT_OUTLIER_WEIGHT,
            Some(&seed),
            &[0.0, 0.20],
        );

        assert!(
            clean
                .iter()
                .zip(&dirty)
                .any(|(a, b)| (a.0 - b.0).abs() > 0.5),
            "the same reads from a clean library and a contaminated one scored alike: \
             {clean:?} against {dirty:?}"
        );
    }
}
