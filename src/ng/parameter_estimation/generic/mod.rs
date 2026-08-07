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
pub mod depth_and_alt_reads;
pub mod depth_bins;
pub mod histogram;
pub mod noise_model;
pub mod read_group_error_rate;
pub mod runs;

use std::collections::BTreeMap;

use smallvec::SmallVec;

use crate::ng::parameter_estimation::fitting::FitTermination;
use crate::ng::parameter_estimation::generic::runs::RunsModelFit;
use crate::ng::parameter_estimation::{Estimate, ParameterEstimationError};
use crate::ng::types::{Bp, ErrorRate, GenotypeFrequency, InbreedingF, Ploidy, ReadGroupId};

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
    /// The genotype frequencies, **one set per ploidy present**. A genome with a
    /// haploid sex chromosome has two entries; today's runs have one.
    pub rates: BTreeMap<Ploidy, Estimate<SampleRates>>,
    /// The fraction of the analysable genome lying in runs of homozygosity.
    /// Diploid-only, and `None` when no diploid region exists: above two copies `F`
    /// needs several identity-by-descent coefficients and is deferred.
    pub inbreeding: Option<Estimate<InbreedingF>>,
    /// What the runs model fitted alongside `F`, when it ran.
    pub runs_model: Option<RunsModelFit>,
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
    pub error_rate: BTreeMap<ReadGroupId, Estimate<ErrorRate>>,
    pub rates: BTreeMap<Ploidy, Estimate<SampleRates>>,
    pub termination: FitTermination,
}

#[cfg(test)]
mod tests {
    use super::*;

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
