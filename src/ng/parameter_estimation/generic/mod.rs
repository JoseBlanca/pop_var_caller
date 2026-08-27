//! The SNP/indel path: two tallies of what a sample's sites looked like, the
//! vocabulary they are keyed on, and the four numbers fitted from them.
//!
//! Two accumulators, differing only in how a site is keyed. The **read-group** one
//! enters a site once per read group that covered it, because an error rate describes
//! the chemistry and two libraries of one sample can genuinely differ. The
//! **windowed** one enters that same site once, at its total depth, because
//! heterozygosity describes the individual — one genome has one heterozygosity
//! however many libraries were used to read it. Neither is derivable from the other
//! once a sample has two read groups
//! (`arch/parameter_prepass_generic.md` §3).
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_generic.md` and its architecture
//! companion. The accumulators and the fits land across Milestones B, C, E and F;
//! the vocabulary below is Milestone A.

pub mod accumulators;
pub mod calibration;
pub mod coupled_fit;
pub mod depth_and_alt_reads;
pub mod depth_bins;
pub mod estimate;
#[cfg(test)]
mod expected_counts;
pub mod fallback;
pub mod histogram;
pub mod noise_model;
pub mod read_group_error_rate;
#[cfg(test)]
mod real_alignments;
#[cfg(test)]
mod recovery;
pub mod runs;
#[cfg(test)]
mod truth_anchors;

use std::collections::{BTreeMap, BTreeSet};

use smallvec::SmallVec;

use crate::ng::parameter_estimation::fitting::FitTermination;
use crate::ng::parameter_estimation::generic::calibration::MintedReadErrors;
use crate::ng::parameter_estimation::generic::runs::RunsModelFit;
use crate::ng::parameter_estimation::{Estimate, ParameterEstimationError};
use crate::ng::types::{
    Bp, DomainError, ErrorRate, GenotypeFrequency, InbreedingF, Ploidy, ReadGroupId,
};

/// How often a site is one where reads disagree with the reference far more than the
/// library's chemistry explains, and how badly they disagree there.
///
/// **Why the generic path needs a second class of site at all.** Its noise model is one
/// substitution rate per read group, and measured on HG002's confident regions the body of
/// that distribution is right while its tail is not: 818 loci carrying no benchmark variant
/// showed three or more alternative reads where the model predicts 29. Mismapped reads and
/// error-prone sequence contexts produce sites like that, and a mixture over three genotypes
/// has exactly one class that can explain them — so the surplus arrived as heterozygosity,
/// **1.41 times the benchmark count on that sample**
/// (`research/noise_model_overdispersion_2026-08-10.md`).
///
/// A site is *clean* with probability `1 − noisy_fraction` and *noisy* with probability
/// `noisy_fraction`; the genotype emission then uses that site's own rate. Fitted
/// independently at two depths on HG002 this comes out at about one site in 110 disagreeing
/// at 4–5% against a clean 0.19%.
///
/// **What makes a site noisy, and why more than one cause matters here** (owner,
/// 2026-08-10). At least three things produce a population of sites where reads disagree with
/// the reference far more often than chemistry explains, and they do not belong to the same
/// thing:
///
/// - **Duplications the reference does not carry.** A sequenced sample holding two copies of
///   a region the reference holds once collects both copies' reads at one locus, and the
///   positions where the copies differ show alternative reads at every depth. This is a
///   property of the **genome**, so it is shared by every library made from it.
/// - **Contamination in a library.** Reads from another individual raise the alternative
///   count at exactly the loci where that individual differs from this one. This is a
///   property of the **library**, and two libraries of the same sample can differ in it.
/// - **Error-prone sequence context and mismapping**, which is partly the library's too:
///   mapping difficulty depends on read length and insert size.
///
/// So a noisy population is expected in most samples rather than in unlucky ones, and its
/// share is not necessarily the same in two libraries of one sample.
///
/// **Fitted per sample and shared across its read groups all the same, while `ε` stays per
/// read group.** No data distinguishes the two: every sample in both cohorts carries one
/// library, so a per-library share would be fitted from the same sites as the per-sample one.
/// The per-sample choice is also the one that keeps each library's rate on a one-dimensional
/// ladder. **An assumption to revisit as soon as a multi-library alignment exists**, not a
/// measured conclusion — and contamination is the case that would break it first, because it
/// can raise one library's noisy share while leaving its sibling's untouched.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct SiteNoise {
    noisy_fraction: f64,
    noisy_error_rate: ErrorRate,
}

impl SiteNoise {
    /// The only constructor.
    ///
    /// # Errors
    ///
    /// [`DomainError::SiteNoiseFraction`] when the fraction is not a probability. The rate
    /// is already checked by its own type.
    pub fn try_new(noisy_fraction: f64, noisy_error_rate: ErrorRate) -> Result<Self, DomainError> {
        if (0.0..=1.0).contains(&noisy_fraction) {
            Ok(Self {
                noisy_fraction,
                noisy_error_rate,
            })
        } else {
            Err(DomainError::SiteNoiseFraction(noisy_fraction))
        }
    }

    /// The share of sites drawn from the noisy class.
    #[must_use]
    pub fn noisy_fraction(self) -> f64 {
        self.noisy_fraction
    }

    /// How often a read disagrees with the reference at a noisy site.
    #[must_use]
    pub fn noisy_error_rate(self) -> ErrorRate {
        self.noisy_error_rate
    }

    /// The rate a read disagrees at a site drawn at **random** — the two classes' rates
    /// weighted by how often each class occurs.
    ///
    /// **This is what a sample emits as its error rate, and the choice is deliberate**
    /// (owner-approved, 2026-08-10). It keeps `Estimate<ErrorRate>` and every consumer of it
    /// unchanged; it is the quantity the model-free count at benchmark
    /// homozygous-reference positions measures, so `arch/parameter_prepass_generic.md` §9's
    /// anchor still applies — measured 2.344 × 10⁻³ against a model-free 2.263 × 10⁻³, 3.6%
    /// high and inside one rung of the error-rate ladder. Emitting the *clean* rate instead
    /// would report 16% **below** the model-free count, which that same section calls an
    /// unambiguous bug.
    ///
    /// # Panics
    ///
    /// Never through this type: a convex combination of two values in `[0, 1]` is in
    /// `[0, 1]`, so the checked constructor cannot reject it.
    #[must_use]
    pub fn marginal_error_rate(self, clean: ErrorRate) -> ErrorRate {
        let marginal = (1.0 - self.noisy_fraction) * clean.get()
            + self.noisy_fraction * self.noisy_error_rate.get();
        ErrorRate::try_new(marginal)
            .expect("a convex combination of two probabilities is a probability")
    }
}

/// Which fixed-width window of the reference a locus falls in — its start position
/// divided by [`INBREEDING_WINDOW_BP`], within a contig. Windows never span contigs.
///
/// Unconstrained: any `u32` is a legal window number, so the field is public and
/// there is no checked constructor — the same call
/// [`ContigId`](crate::ng::types::ContigId) makes.
///
/// It stays in this module, and not in the shared vocabulary, because the window
/// exists to serve the runs model that fits the inbreeding coefficient and nothing
/// else.
///
/// **A note for whoever writes the division** (Milestone C): ng's
/// [`Position`](crate::ng::types::Position) is **1-based**, so a naive
/// `start / INBREEDING_WINDOW_BP` puts positions 1–99,999 in window 0 and gives every
/// later window 100,000 bases. Deciding what to do about the first window's 99,999 is
/// that step's, not this type's.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct WindowIndex(pub u32);

impl WindowIndex {
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

/// The width of an inbreeding window.
///
/// **Fixed, not a knob**: a window size is not a quantity a user is in a position to
/// choose, and an unsettable knob is worse than a constant because it invites a wrong
/// answer and offers no way to recognise one
/// (`spec/parameter_prepass_generic.md` §4).
///
/// It is the grain the runs model classifies as inside or outside a run of
/// homozygosity, and 100 kb is set by the shortest run worth resolving — about
/// 300 kb — not by what the accumulator costs. Both organisms in view are far above
/// that: a tomato landrace is homozygous over tens of megabases, and a consanguineous
/// human's segments run from 5 to 50 Mb.
pub const INBREEDING_WINDOW_BP: Bp = Bp(100_000);

// ---------------------------------------------------------------------
// The error-rate ladder — the candidate rates the **profile scan** steps
// through. That scan scores every candidate error rate in turn, refitting
// the genotype frequencies at each and keeping the best-scoring rung: a
// profile likelihood (`arch/parameter_prepass_generic.md` §4.2).
//
// **Phred appears in step 4 only here.** The rungs themselves are
// probabilities; the Phred scale is how the ladder is *spaced*, because
// these rates span orders of magnitude and the distance that matters
// between two of them is a ratio rather than a difference.
//
// There is deliberately no newtype for a Phred-scaled *rate*. `types.rs`
// carries Phred only as the integer read qualities `BaseQual` and
// `MapQual`, and it already has `LogProb` for the logarithm of a
// probability; a second log-scaled probability type in a different base
// would make a base mix-up a plausible wrong number instead of a compile
// error — the very hazard `LogProb` exists to prevent
// (`arch/parameter_prepass_generic.md` §2.1).
// ---------------------------------------------------------------------

/// The ladder's **noisiest** rung: Phred 10, an error rate of 0.1 — one base in ten
/// wrong, which is worse than any usable run.
///
/// **Fixed, not a knob**, and the range is DRAGstr's own for the same kind of grid
/// (`spec/parameter_prepass.md` §3). A read group whose true rate lies outside
/// Phred 10–50 — a bad run, heavy contamination — has its answer clamped to an edge,
/// and the remedy is not a wider ladder: the scan reports an endpoint argmax so the
/// railed fit announces itself (`arch/parameter_prepass_generic.md` §4.2).
pub const ERROR_RATE_LADDER_MIN_PHRED: f32 = 10.0;

/// The ladder's **cleanest** rung: Phred 50, an error rate of 0.00001 — one base in a
/// hundred thousand, below what any current chemistry delivers.
///
/// **Fixed, not a knob**, for the reason [`ERROR_RATE_LADDER_MIN_PHRED`] gives; the
/// two edges are one decision.
pub const ERROR_RATE_LADDER_MAX_PHRED: f32 = 50.0;

/// The spacing between rungs: a quarter of a Phred, so adjacent rungs differ by a
/// factor of `10^0.025` — about 6% — in probability.
///
/// The spec argues that is below what a caller can feel, and **marks the argument
/// soft**: "a few percent" is an argument from what a prior does, not a measurement
/// of what a caller tolerates, and it is untested until the synthetic fits run
/// (`spec/parameter_prepass.md` §3). On that argument the scan is a single flat pass
/// with no refinement stage.
pub const ERROR_RATE_LADDER_STEP_PHRED: f32 = 0.25;

/// How many rungs the ladder has.
///
/// **Stated, not derived.** An earlier version computed the count by rounding
/// `(max − min) / step` and casting to an integer, which fails silently two ways: a
/// maximum off the step grid is absorbed by the rounding, so the ladder stops short
/// of the constant named `MAX`; and an inverted pair gives a negative float whose
/// `as u32` cast **saturates to zero**, collapsing the ladder to a single rung. A
/// one-rung ladder would set the endpoint-argmax flag for every read group, which is
/// the one bit standing between a railed fit and a plausible-looking number.
///
/// Stating it here and checking it against the three Phred constants — at build time
/// for the ordering, by test for the arithmetic — leaves neither failure silent.
pub const ERROR_RATE_LADDER_RUNGS: usize = 161;

// The ladder runs upward from a non-negative Phred in positive steps. Checked at
// build time rather than by test, because these three are `pub const`s a later edit
// can change, and `error_rate_ladder`'s `PANIC-FREE` claim below rests on them.
const _: () = assert!(
    ERROR_RATE_LADDER_MIN_PHRED >= 0.0
        && ERROR_RATE_LADDER_MAX_PHRED > ERROR_RATE_LADDER_MIN_PHRED
        && ERROR_RATE_LADDER_STEP_PHRED > 0.0,
    "the error-rate ladder runs upward from a non-negative Phred in positive steps"
);

/// The error rates the profile scan steps through: [`ERROR_RATE_LADDER_RUNGS`] rungs
/// from [`ERROR_RATE_LADDER_MIN_PHRED`] upward in steps of
/// [`ERROR_RATE_LADDER_STEP_PHRED`] — **ascending in Phred**, and so descending in
/// probability, from 0.1 down to 0.00001.
///
/// Built rather than stored as a table, so the constants above are the single
/// statement of the ladder's shape. It allocates a fresh vector each call; the scan
/// builds it once per fit and re-walks the slice at every rung, so nothing calls this
/// in a loop.
#[must_use]
pub fn error_rate_ladder() -> Vec<ErrorRate> {
    let min_phred = f64::from(ERROR_RATE_LADDER_MIN_PHRED);
    let step_phred = f64::from(ERROR_RATE_LADDER_STEP_PHRED);
    (0..ERROR_RATE_LADDER_RUNGS)
        .map(|rung| {
            let phred = min_phred + rung as f64 * step_phred;
            // PANIC-FREE: the const assertion above pins the ladder to non-negative
            // Phred values ascending in positive steps, so `10^(-phred/10)` is in
            // `(0, 1]` and the checked constructor cannot reject it.
            ErrorRate::try_new(10f64.powf(-phred / 10.0)).unwrap_or_else(|rejected| {
                panic!("ladder rung {rung} (Phred {phred}) is not a probability: {rejected}")
            })
        })
        .collect()
}

// ---------------------------------------------------------------------
// When there is not enough data — the floors, and the one default.
// ---------------------------------------------------------------------

/// Fewest sites a fit will accept before it borrows or fails.
///
/// **Soft** — a guard against fitting a rate from a handful of sites and emitting it as
/// though it were measured, not a precision target. The spec's own precision figure is
/// about a different regime and does not bear on this number: six million read
/// observations pin an error rate near 0.001 to about one part in eighty, and at three
/// reads a site that is two million sites, 200 times this floor
/// (`spec/parameter_prepass.md` §4.1). What a fit at 10,000 sites is worth was not
/// measured.
pub const MIN_SITES_TO_FIT: u64 = 10_000;

/// The per-base error rate used when none can be fitted and none was supplied.
///
/// **Soft, and the only defaulted parameter on this path.** Chemistry varies far less
/// between runs than biology does between samples, so a stated constant is defensible
/// here in a way it is not for heterozygosity or for inbreeding — which is why those
/// two fail instead (`arch/parameter_prepass_generic.md` §5.4).
pub const DEFAULT_ERROR_RATE: f64 = 0.001;

/// The cap on the alternation between the error rates and the genotype frequencies.
///
/// **Soft, and generous against what the loop needs.** The alternation converges
/// linearly and was measured reaching the truth in all 25 worlds tried from a start at
/// three times the true rates and half the true frequencies
/// (`spec/parameter_prepass_generic.md` §5.1); a single-library sample takes **one**
/// iteration, because the two tables it reads are then the same table. The cap exists
/// because the outer alternation — unlike the concave inner climb — has no convergence
/// proof, so it is capped, the best-scoring iterate is kept, and the termination is
/// reported rather than a stalled fit arriving silently.
pub const MAX_COUPLED_FIT_ITERATIONS: u32 = 20;

// ---------------------------------------------------------------------
// What the generic path emits.
// ---------------------------------------------------------------------

/// A sample's genotype frequencies at one ploidy: one entry per number of alternative
/// copies, `0..=P`, summing to one.
///
/// A vector rather than two named fields, because at `P = 4` there are five entries and
/// the intermediate dosages have no diploid name. The two a diploid caller reads are
/// the accessors below.
///
/// **The fields are private and the constructor is checked, because the accessors read
/// by dosage.** `homozygous_non_reference_rate()` returns the *last* entry, so a
/// ploidy-2 set holding one entry would hand back the homozygous-**reference** rate —
/// near 1.0 — under the homozygous-**non-reference** name, where the truth is near
/// 0.001. That is a wrong number with no symptom, in the module whose whole difficulty
/// is that wrong numbers here have none. Guarding the accessors instead was the other
/// option and it is worse: it answers `None` for a malformed set, which no caller can
/// act on and which the doc below defines as meaning something else entirely.
#[derive(Clone, PartialEq, Debug)]
pub struct SampleRates {
    ploidy: Ploidy,
    by_alt_copies: SmallVec<[GenotypeFrequency; 5]>,
}

impl SampleRates {
    /// The only constructor. `by_alt_copies` must hold one frequency per dosage
    /// `0..=ploidy`, summing to one.
    ///
    /// **Fallible, though the caller is our own fit.** A set off the simplex means the
    /// climb that produced it is broken and there is nothing a caller could do — so the
    /// fits construct through this door and `.expect()`, the same shape the four
    /// constrained scalars use. What makes it worth a check rather than a comment is
    /// that the failure is silent: see the type's doc.
    ///
    /// The sum is compared with a tolerance, because the frequencies arrive from a
    /// floating-point climb over the simplex and an exact `== 1.0` would reject a
    /// correct fit.
    pub fn try_new(
        ploidy: Ploidy,
        by_alt_copies: SmallVec<[GenotypeFrequency; 5]>,
    ) -> Result<Self, ParameterEstimationError> {
        let total: f64 = by_alt_copies.iter().map(|f| f.get()).sum();
        let one_per_dosage = by_alt_copies.len() == usize::from(ploidy.get()) + 1;

        if !one_per_dosage || (total - 1.0).abs() > SIMPLEX_TOLERANCE {
            return Err(ParameterEstimationError::GenotypeFrequenciesOffSimplex {
                ploidy,
                entries: by_alt_copies.len(),
                total,
            });
        }
        Ok(Self {
            ploidy,
            by_alt_copies,
        })
    }

    /// How many copies of the genome this set describes.
    #[must_use]
    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    /// Every frequency, indexed by how many of the individual's copies are
    /// non-reference: entry 0 is the homozygous-reference rate, the last entry the
    /// homozygous-non-reference one.
    #[must_use]
    pub fn by_alt_copies(&self) -> &[GenotypeFrequency] {
        &self.by_alt_copies
    }

    /// How often the individual's two copies of a site differ — `Hobs` in the
    /// population genetics.
    ///
    /// **Diploid only.** Above two copies a site can carry one, two or three
    /// alternative copies, and calling all of them "heterozygous" throws the dosage
    /// away; the replacement is gene diversity — the chance that two copies drawn at
    /// random from the individual differ — which reduces to this at `P = 2` and is
    /// deferred. `None` says "this genome does not have one of these", not "it was not
    /// measured", and the checked constructor is what keeps that the only reading.
    #[must_use]
    pub fn observed_heterozygosity(&self) -> Option<GenotypeFrequency> {
        (self.ploidy.get() == 2).then(|| self.by_alt_copies[1])
    }

    /// How often **every** copy is non-reference: the last entry, at any ploidy.
    ///
    /// **Not a leftover of heterozygosity — it measures something else.** How often an
    /// individual carries a non-reference allele at all is *heterozygosity plus this
    /// rate*, and it is that **sum** that belongs to the pair (individual, reference):
    /// swap in a different accession as the reference and it changes, while
    /// heterozygosity and inbreeding do not. The two terms also come apart in the
    /// direction that matters here — a selfing landrace far from the reference is mostly
    /// homozygous and mostly non-reference at once, so a caller whose prior assumes
    /// "non-reference implies rare" is wrong on exactly that sample, and tomato's
    /// reference is one cultivated accession.
    #[must_use]
    pub fn homozygous_non_reference_rate(&self) -> GenotypeFrequency {
        // PANIC-FREE: `try_new` is the only constructor and it rejects a set with fewer
        // than `ploidy + 1` entries, so there is always a last one.
        self.by_alt_copies[self.by_alt_copies.len() - 1]
    }
}

/// How far the sum of a genotype-frequency set may sit from one. The frequencies arrive
/// from a floating-point climb over the simplex, so an exact comparison would reject a
/// correct fit; this is loose enough for accumulated rounding over five entries and far
/// too tight to admit a set that is genuinely off.
const SIMPLEX_TOLERANCE: f64 = 1e-9;

/// Everything the generic path estimates for one sample.
///
/// Named for what it holds rather than "priors": half of it is noise-model terms, and
/// what a caller builds into a prior is the caller's design.
#[derive(Clone, PartialEq, Debug)]
pub struct GenericSampleParameters {
    /// One per read group. Chemistry belongs to the library rather than to the
    /// individual and does not vary with ploidy, so there is one of these however many
    /// ploidies the genome holds.
    pub error_rate: BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    /// **How wrong this sample's reads said they were, summed per read group** — the
    /// denominator the calling step divides [`Self::error_rate`] by
    /// (`doc/devel/ng/spec/read_likelihoods.md` §3.2).
    ///
    /// A read arrives carrying its own claim about how likely it is to be wrong, from its base
    /// and mapping qualities. The fit above says how often a read of this library *actually*
    /// disagrees. The calling step scores each read at its own claim multiplied by the ratio
    /// between the two, so it needs both halves — and this is the half that is a sum over the
    /// reads rather than a fit, which is why it comes off the tally
    /// ([`GenericAccumulators::minted_errors`](accumulators::GenericAccumulators::minted_errors))
    /// unchanged rather than out of a fit.
    ///
    /// **Carried here because the tally does not outlive the fit.** The streaming entry point
    /// drops its accumulator as soon as the fit returns, so this is the last moment the totals
    /// can be reached — and until 2026-08-27 they went with it, leaving nothing that assembled a
    /// run's calling parameters able to supply the denominator at all. That much does not fail
    /// quietly: `RunParameters::assemble` refuses a fitted rate whose read group has no total
    /// here, and the reverse, so a run built without them stops at assembly naming the first read
    /// group. **What is quiet is a total under the wrong key** — a map with the right identifiers
    /// and another read group's numbers in them, which no assertion can see and which moves every
    /// scale it touches.
    ///
    /// **A read group with no entry is one that put no complete observation anywhere**, which
    /// is not the same as one whose reads all read perfectly: a sum of zero over reads that
    /// were seen is an ordinary value, since a read at Phred 0 contributes `ln 1 = 0`. The
    /// calling step's own assembly refuses a fitted rate whose read group has no total here,
    /// and the reverse, because the two are one pass over one set of reads.
    pub minted_errors: BTreeMap<ReadGroupId, MintedReadErrors>,
    /// The genotype frequencies, **one set per ploidy present**. A genome with a
    /// haploid sex chromosome has two entries; today's runs have one.
    pub rates: BTreeMap<Ploidy, Estimate<SampleRates>>,
    /// The fraction of the analysable genome lying in runs of homozygosity.
    /// Diploid-only, and `None` when no diploid region exists: above two copies `F`
    /// needs several identity-by-descent coefficients and is deferred.
    pub inbreeding: Option<Estimate<InbreedingF>>,
    /// What the runs model fitted alongside `F`, when it ran.
    pub runs_model: Option<RunsModelFit>,
    /// The sample's second class of site, when its data asked for one — as
    /// [`CoupledFit::site_noise`], and already folded into `error_rate`.
    pub site_noise: Option<SiteNoise>,
    /// **This sample wanted a second class of site outside the model's range and was refused**
    /// — see [`CoupledFit::site_noise_off_the_ladder`]. `site_noise` is then `None` for a
    /// reason a consumer should not read as *no second class was needed*.
    pub site_noise_off_the_ladder: bool,
    /// **The read groups whose error rate was clamped to an end of the ladder** — as
    /// [`CoupledFit::error_rate_on_a_ladder_end`]. Empty on every sample measured so far.
    pub error_rate_on_a_ladder_end: BTreeSet<ReadGroupId>,
    /// How the coupled error-rate/frequency fit ended.
    pub coupled_fit: FitTermination,
}

/// What the coupled fit returns: the two quantities that had to be fitted together, and
/// how the alternation between them ended.
///
/// **One call for the whole sample, not one per ploidy.** The error rates come out per
/// read group and span every ploidy that group covered — chemistry does not know about
/// chromosomes — while the genotype frequencies come out one set per ploidy present,
/// because a haploid region has two genotype classes and a diploid three.
#[derive(Clone, PartialEq, Debug)]
pub struct CoupledFit {
    /// **The rate a read disagrees at a *clean* site**, per read group — not the rate at a
    /// site drawn at random, which is what a sample emits.
    ///
    /// The two differ by 15% on a sample like HG002, and the distinction is the whole reason
    /// the pair travels: everything inside the fit scores with this rate and [`site_noise`]
    /// beside it, because that is the pair the scoring rule takes, and
    /// [`GenericSampleParameters::error_rate`] is where the two are folded into the one
    /// number a consumer reads. Folding sooner would put the tail misspecification back
    /// inside the runs model.
    ///
    /// Where the fit found no second class the two are the same number.
    ///
    /// [`site_noise`]: CoupledFit::site_noise
    pub error_rate: BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    /// **The read groups whose rate was clamped to an end of the ladder rather than found
    /// inside it** — Phred 10 or Phred 50, one base in ten or one in a hundred thousand.
    ///
    /// An argmax on an edge is the edge of the search and not a maximum in it, which
    /// `arch/parameter_prepass_generic.md` §9 names as one of the two ways this estimator
    /// returns a confident wrong number rather than failing. **It was computed and then
    /// dropped**, so no consumer could tell a clamped rate from a fitted one; the second class
    /// of site reports the same shape in [`site_noise_off_the_ladder`], and this is the other
    /// half of it.
    ///
    /// Empty on every sample measured so far: a real library's rate sits far inside Phred 10
    /// to 50. Only groups whose rate was **fitted here** can appear — a borrowed, supplied or
    /// defaulted rate is not the argmax of anything.
    ///
    /// [`site_noise_off_the_ladder`]: CoupledFit::site_noise_off_the_ladder
    pub error_rate_on_a_ladder_end: BTreeSet<ReadGroupId>,
    pub rates: BTreeMap<Ploidy, Estimate<SampleRates>>,
    /// The sample's second class of site, when its data asked for one.
    ///
    /// **The other half of a pair, not a diagnostic to be ignored here.** Unlike
    /// [`GenericSampleParameters::site_noise`], whose partner rate has already been
    /// marginalised, this one is still needed to say what [`error_rate`] means: the two
    /// together are the sample's noise model, and reading the rate without the share
    /// understates how often a read disagrees.
    ///
    /// [`error_rate`]: CoupledFit::error_rate
    pub site_noise: Option<SiteNoise>,
    /// **The sample asked for a second class of site this model does not cover, and was
    /// refused** — so [`site_noise`] is `None` for a reason quite unlike the ordinary one.
    ///
    /// The error-rate ladder runs Phred 10 to 50 because that is the range of *sequencing*
    /// noise (`spec/parameter_prepass.md` §3). A best second class sitting on its coarsest
    /// rung is the search asking to leave that range, and what such a sample holds is a
    /// population of positions the model does not describe — duplications the reference does
    /// not carry look like this, since about half the reads disagree where the copies differ,
    /// five times that rung. Owner's call, 2026-08-10: refuse it and fit one rate, rather than
    /// widen the model until it stops serving the samples that do meet its assumptions.
    ///
    /// **Without this bit the two cases are indistinguishable**, and they are not alike: an
    /// ordinary `None` says the data did not ask for a second class, and this one says it
    /// asked for something out of range. On the five real alignments run so far it is set on
    /// two — tomato SRR7279482 and SRR7279483.
    ///
    /// [`site_noise`]: CoupledFit::site_noise
    pub site_noise_off_the_ladder: bool,
    pub termination: FitTermination,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rate(x: f64) -> ErrorRate {
        ErrorRate::try_new(x).expect("a probability")
    }

    /// The share of noisy sites is a probability and nothing else is accepted, including
    /// the two values a floating-point fit can drift to.
    #[test]
    fn a_noisy_site_share_outside_zero_to_one_is_refused() {
        assert!(
            SiteNoise::try_new(0.0, rate(0.05)).is_ok(),
            "no noisy sites at all"
        );
        assert!(
            SiteNoise::try_new(1.0, rate(0.05)).is_ok(),
            "every site noisy"
        );
        for refused in [
            -1e-9,
            1.000_000_001,
            f64::NAN,
            f64::INFINITY,
            f64::NEG_INFINITY,
        ] {
            assert!(
                matches!(
                    SiteNoise::try_new(refused, rate(0.05)),
                    Err(DomainError::SiteNoiseFraction(_))
                ),
                "{refused} is not a share of sites"
            );
        }
    }

    /// **The number a sample emits, worked by hand.** At the fitted HG002 30x values —
    /// 0.88% of sites noisy at 5.29e-2, the rest clean at 1.895e-3 — the marginal is
    /// 0.9912 x 1.895e-3 + 0.0088 x 5.29e-2 = 1.878324e-3 + 4.6552e-4 = 2.343844e-3, which
    /// is what the model-free count
    /// of 2.263e-3 is compared against (3.6% high, inside one ladder rung). Pinned here
    /// because a transposition of the two rates, or of the share and its complement, still
    /// returns a probability and would be caught by nothing else.
    #[test]
    fn the_emitted_rate_is_the_two_classes_weighted_by_how_often_each_occurs() {
        let measured = SiteNoise::try_new(0.0088, rate(5.29e-2)).expect("a share and a rate");
        let marginal = measured.marginal_error_rate(rate(1.895e-3)).get();
        assert!(
            (marginal - 2.343_844e-3).abs() < 1e-9,
            "expected 2.343844e-3 by hand, got {marginal:e}"
        );

        // Transposing the two rates gives a wildly different answer, so the test bites.
        let transposed = SiteNoise::try_new(0.0088, rate(1.895e-3))
            .expect("a share and a rate")
            .marginal_error_rate(rate(5.29e-2))
            .get();
        assert!(
            transposed > 20.0 * marginal,
            "the transposition must not be mistakable for the truth: {transposed:e}"
        );
    }

    /// **No noisy sites means the emitted rate is the clean one, exactly** — the property
    /// that lets a one-class sample keep every number it has today, and the reason the
    /// scoring rule can branch on it.
    #[test]
    fn a_sample_with_no_noisy_sites_emits_its_clean_rate_unchanged() {
        let none = SiteNoise::try_new(0.0, rate(0.5)).expect("a share and a rate");
        for clean in [1e-5, 1.895e-3, 0.1] {
            assert_eq!(
                none.marginal_error_rate(rate(clean)).get(),
                clean,
                "a sample with no noisy sites reports its clean rate"
            );
        }
    }

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("a positive copy number")
    }

    fn rates(copies: u8, by_alt_copies: &[f64]) -> SampleRates {
        SampleRates::try_new(ploidy(copies), frequencies(by_alt_copies))
            .expect("a well-formed set for this ploidy")
    }

    fn frequencies(values: &[f64]) -> SmallVec<[GenotypeFrequency; 5]> {
        values
            .iter()
            .map(|&f| GenotypeFrequency::try_new(f).expect("a frequency in [0, 1]"))
            .collect()
    }

    /// A diploid sample reads its heterozygosity off the middle entry, and the three
    /// frequencies sum to one — the simplex the fit climbs on.
    #[test]
    fn a_diploid_sample_reads_its_two_rates_off_three_frequencies_summing_to_one() {
        let diploid = rates(2, &[0.9885, 0.0105, 0.0010]);

        let total: f64 = diploid.by_alt_copies().iter().map(|f| f.get()).sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "the frequencies sum to one: {total}"
        );

        assert_eq!(
            diploid.observed_heterozygosity().map(|h| h.get()),
            Some(0.0105)
        );
        assert_eq!(diploid.homozygous_non_reference_rate().get(), 0.0010);
    }

    /// **Above two copies there is no heterozygosity to read**, and the accessor says so
    /// rather than inventing one. A tetraploid site can carry one, two or three
    /// alternative copies; calling all three "heterozygous" would throw the dosage away,
    /// and the replacement — gene diversity — is deferred. The homozygous-non-reference
    /// rate is defined at every ploidy, because "every copy differs" always means one
    /// thing.
    #[test]
    fn heterozygosity_is_absent_above_diploidy_where_the_homozygous_rate_is_not() {
        let tetraploid = rates(4, &[0.970, 0.018, 0.007, 0.003, 0.002]);

        let total: f64 = tetraploid.by_alt_copies().iter().map(|f| f.get()).sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "the frequencies sum to one: {total}"
        );

        assert_eq!(tetraploid.by_alt_copies().len(), 5, "one entry per 0..=P");
        assert_eq!(tetraploid.observed_heterozygosity(), None);
        assert_eq!(tetraploid.homozygous_non_reference_rate().get(), 0.002);
    }

    /// A haploid region has two genotype classes, not three, and it is not diploid — so
    /// it has no heterozygosity either.
    ///
    /// **The boundary is tested on both sides because the accessor's condition is an
    /// equality**, and the two ways to loosen it fail in opposite directions: `>= 2`
    /// lets a *tetraploid* answer, which the tetraploid test above catches, and `<= 2`
    /// lets a haploid answer, which only this one does.
    #[test]
    fn a_haploid_region_has_two_classes_and_no_heterozygosity() {
        let haploid = rates(1, &[0.9985, 0.0015]);

        assert_eq!(haploid.by_alt_copies().len(), 2);
        assert_eq!(haploid.observed_heterozygosity(), None);
        assert_eq!(haploid.homozygous_non_reference_rate().get(), 0.0015);
    }

    /// **A set with the wrong number of entries for its ploidy is rejected, and this is
    /// the test that matters most in this file.** The accessors read *by dosage*: before
    /// the constructor was checked, a ploidy-2 set holding only entry 0 answered
    /// `homozygous_non_reference_rate()` with that entry — the homozygous-*reference*
    /// rate, near 1.0, under the homozygous-*non-reference* name, where the truth is
    /// near 0.001.
    #[test]
    fn a_set_with_the_wrong_number_of_entries_for_its_ploidy_is_rejected() {
        for (copies, values) in [
            (2u8, &[1.0][..]),                  // the inversion above
            (2, &[0.5, 0.3, 0.1, 0.1][..]),     // one entry too many
            (4, &[0.5, 0.5][..]),               // a tetraploid keyed as a haploid
            (1, &[0.9885, 0.0105, 0.0010][..]), // a haploid keyed as a diploid
        ] {
            let rejected = SampleRates::try_new(ploidy(copies), frequencies(values));
            assert!(
                matches!(
                    rejected,
                    Err(ParameterEstimationError::GenotypeFrequenciesOffSimplex { .. })
                ),
                "ploidy {copies} with {} entries should not construct",
                values.len()
            );
        }
    }

    /// A set of the right length whose entries do not sum to one is rejected too: that
    /// is not a distribution, and every fit that reads one assumes it is. The tolerance
    /// admits a floating-point climb's accumulated rounding and nothing wider.
    #[test]
    fn a_set_that_does_not_sum_to_one_is_rejected_within_a_rounding_tolerance() {
        assert!(SampleRates::try_new(ploidy(2), frequencies(&[0.5, 0.3, 0.1])).is_err());
        assert!(SampleRates::try_new(ploidy(2), frequencies(&[0.5, 0.4, 0.4])).is_err());

        // A third of the genome in each class: the entries sum to one only after
        // rounding, and a fit is entitled to produce this.
        let third = 1.0 / 3.0;
        assert!(SampleRates::try_new(ploidy(2), frequencies(&[third, third, third])).is_ok());
    }

    /// The floors are stated as tests because each is a number a later reader will want
    /// to change, and each has a measurement or a stated decision behind it rather than
    /// an argument.
    #[test]
    fn the_fit_floors_are_the_measured_ones() {
        assert_eq!(MIN_SITES_TO_FIT, 10_000);
        assert_eq!(DEFAULT_ERROR_RATE, 0.001);
        assert_eq!(MAX_COUPLED_FIT_ITERATIONS, 20);
    }

    /// The default error rate has to be a rate the type accepts, or the fallback rung
    /// that exists for a sample with no fittable read group cannot be taken at all.
    #[test]
    fn the_default_error_rate_is_a_constructible_rate() {
        assert_eq!(
            ErrorRate::try_new(DEFAULT_ERROR_RATE).map(|r| r.get()),
            Ok(0.001)
        );
    }

    /// The ladder's shape is a correctness parameter of the error-rate fit, not a
    /// convenience: 161 rungs at quarter-Phred spacing is what the spec argues is
    /// finer than a caller can feel, and the whole "no refinement stage" decision
    /// rests on it.
    #[test]
    fn the_error_rate_ladder_spans_phred_10_to_50_in_161_rungs() {
        let ladder = error_rate_ladder();

        assert_eq!(ladder.len(), 161, "(50 - 10) / 0.25 + 1");
        assert!(
            (ladder[0].get() - 0.1).abs() < 1e-12,
            "first rung is Phred 10, got {}",
            ladder[0].get()
        );
        assert!(
            (ladder[160].get() - 1e-5).abs() < 1e-17,
            "last rung is Phred 50, got {}",
            ladder[160].get()
        );
    }

    /// The rungs are **derived from** the three Phred constants, not merely equal to
    /// them today. Without this, a step that does not divide the span leaves the top
    /// rung short of `ERROR_RATE_LADDER_MAX_PHRED` — so the scan's finest candidate
    /// rate is wrong, and every "railed at the ladder's end" flag downstream is
    /// measured against the wrong edge, with nothing in the output to show it.
    #[test]
    fn the_error_rate_ladder_ends_at_the_phred_constants_it_is_built_from() {
        let ladder = error_rate_ladder();
        let noisiest = 10f64.powf(-f64::from(ERROR_RATE_LADDER_MIN_PHRED) / 10.0);
        let cleanest = 10f64.powf(-f64::from(ERROR_RATE_LADDER_MAX_PHRED) / 10.0);

        assert_eq!(ladder.len(), ERROR_RATE_LADDER_RUNGS);
        assert!((ladder[0].get() - noisiest).abs() <= noisiest * 1e-12);
        let last = ladder.last().expect("the ladder is never empty").get();
        assert!(
            (last - cleanest).abs() <= cleanest * 1e-12,
            "last rung {last} vs {cleanest}"
        );
    }

    /// `ERROR_RATE_LADDER_RUNGS` is stated rather than computed, so something has to
    /// check it against the constants it claims to summarise. A step that leaves a
    /// fractional rung count is the case that would otherwise pass unnoticed.
    #[test]
    fn the_ladder_constants_divide_into_a_whole_number_of_rungs() {
        let steps = f64::from(ERROR_RATE_LADDER_MAX_PHRED - ERROR_RATE_LADDER_MIN_PHRED)
            / f64::from(ERROR_RATE_LADDER_STEP_PHRED);
        assert!(
            (steps - steps.round()).abs() < 1e-6,
            "the ladder's span must be a whole number of steps, got {steps}"
        );
        assert_eq!(steps.round() as usize + 1, ERROR_RATE_LADDER_RUNGS);
    }

    /// Ascending in Phred means **descending** in probability, and the ratio between
    /// neighbours is constant — that is what makes "one rung" a meaningful unit of
    /// distance for the coupled fit, which reports its movement in rungs.
    #[test]
    fn the_error_rate_ladder_rungs_are_a_constant_ratio_apart() {
        let ladder = error_rate_ladder();
        let expected_ratio = 10f64.powf(f64::from(ERROR_RATE_LADDER_STEP_PHRED) / 10.0);

        assert_eq!(
            ladder.windows(2).count(),
            ERROR_RATE_LADDER_RUNGS - 1,
            "a ladder short enough to make the loop below vacuous is itself the bug"
        );
        for pair in ladder.windows(2) {
            let (higher_rate, lower_rate) = (pair[0].get(), pair[1].get());
            assert!(
                lower_rate < higher_rate,
                "rates descend as Phred ascends: {lower_rate} should be below {higher_rate}"
            );
            let ratio = higher_rate / lower_rate;
            assert!(
                (ratio - expected_ratio).abs() < 1e-12,
                "adjacent rungs differ by 10^0.025, got {ratio}"
            );
        }
    }

    /// A window is 100 kb. Stated as a test because the number is load-bearing: the
    /// runs model's noise floor is a function of how many windows a genome has, so
    /// changing this changes the resolution of every reported inbreeding coefficient.
    #[test]
    fn the_inbreeding_window_is_a_hundred_kilobases() {
        assert_eq!(INBREEDING_WINDOW_BP.get(), 100_000);
    }

    #[test]
    fn a_window_index_exposes_its_number() {
        assert_eq!(WindowIndex(7).get(), 7);
    }
}
