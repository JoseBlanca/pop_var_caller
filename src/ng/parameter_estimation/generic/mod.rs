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

pub mod depth_and_alt_reads;
pub mod histogram;
pub mod runs;

use crate::ng::types::{Bp, ErrorRate};

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

#[cfg(test)]
mod tests {
    use super::*;

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
