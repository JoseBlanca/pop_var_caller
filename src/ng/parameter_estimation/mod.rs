//! ng step 4 — the parameters the caller runs on, measured from the sample's own
//! loci before anything is called.
//!
//! Four numbers per sample come out of the SNP/indel path: a per-read-group error
//! rate, the sample's heterozygosity, its homozygous-non-reference rate, and its
//! inbreeding coefficient. They are measured from **every** covered position,
//! including the overwhelming majority that show no alternative allele at all —
//! which is what separates this step from production's estimator, whose caller
//! skips a pure-reference column and so never sees the strongest evidence there is
//! about the error rate (`spec/parameter_prepass.md` §2.1).
//!
//! Design: [`spec/parameter_prepass_generic.md`] (the design and its rationale),
//! [`spec/parameter_prepass.md`] (the shared framing), and
//! [`arch/parameter_prepass_generic.md`] (types and interfaces), all under
//! `doc/devel/ng/`.
//!
//! Two sub-units, split so that the shaping of data and the mathematics on it never
//! live in one file:
//!
//! - [`fitting`] — the mathematics. Knows nothing about markers, loci or windows: it
//!   is given a table of numbers and returns the values that best explain them. A
//!   folder rather than a file because it is the one genuine swappable seam — one
//!   trait, an implementation on this path and a second on the STR path.
//! - [`generic`] — the SNP/indel path: the two accumulators, the cell table, and
//!   what each of the four numbers is fitted from.

pub mod fitting;
pub mod generic;

use crate::ng::types::{Bp, ErrorRate};

/// Which fixed-width window of the reference a locus falls in — `region.start /
/// INBREEDING_WINDOW_BP` within a contig. Windows never span contigs.
///
/// Local to this step and unconstrained, like [`ContigId`](crate::ng::types::ContigId):
/// the window exists to serve the runs model that fits the inbreeding coefficient and
/// nothing else, which is why it is not in the shared vocabulary.
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
/// choose (`spec/parameter_prepass_generic.md` §4). It is the grain the runs model
/// classifies as inside or outside a run of homozygosity, and it is chosen against
/// the runs a genome actually carries, not against a run's data.
pub const INBREEDING_WINDOW_BP: Bp = Bp(100_000);

// ---------------------------------------------------------------------
// The error-rate ladder — the candidate rates the profile scan steps through.
//
// **Phred appears here and nowhere else.** The rungs themselves are
// probabilities; the Phred scale is how the ladder is *spaced*, because
// these rates span orders of magnitude and the distance that matters
// between two of them is a ratio rather than a difference. There is
// deliberately no `Phred` newtype: `types.rs` already has `LogProb` for the
// logarithm of a probability, and a second log-scaled type in a different
// base would make a base mix-up a plausible wrong number instead of a
// compile error — the very hazard `LogProb` exists to prevent
// (`arch/parameter_prepass_generic.md` §2.1).
// ---------------------------------------------------------------------

/// The ladder's coarsest rung: Phred 10, an error rate of 0.1.
pub const ERROR_RATE_LADDER_MIN_PHRED: f32 = 10.0;

/// The ladder's finest rung: Phred 50, an error rate of 0.00001.
pub const ERROR_RATE_LADDER_MAX_PHRED: f32 = 50.0;

/// The spacing between rungs: a quarter of a Phred, so adjacent rungs differ by a
/// factor of `10^0.025` — about 6% — in probability. Finer than a caller can feel
/// (`spec/parameter_prepass.md` §3), which is why the scan is a single flat pass
/// with no refinement stage.
pub const ERROR_RATE_LADDER_STEP_PHRED: f32 = 0.25;

/// The error rates the profile scan steps through: 161 rungs from Phred 10 to
/// Phred 50 in quarter-Phred steps, **ascending in Phred** and so descending in
/// probability, from 0.1 down to 0.00001.
///
/// Built rather than stored as a table so the three constants above are the single
/// statement of the ladder's shape.
pub fn error_rate_ladder() -> Vec<ErrorRate> {
    let span_phred = ERROR_RATE_LADDER_MAX_PHRED - ERROR_RATE_LADDER_MIN_PHRED;
    let top_rung = (span_phred / ERROR_RATE_LADDER_STEP_PHRED).round() as u32;
    (0..=top_rung)
        .map(|rung| {
            let phred = f64::from(ERROR_RATE_LADDER_MIN_PHRED)
                + f64::from(rung) * f64::from(ERROR_RATE_LADDER_STEP_PHRED);
            // PANIC-FREE: `10^(-phred/10)` over a Phred range of 10..=50 is between
            // 0.00001 and 0.1, so it is a finite probability in `[0, 1]` and the
            // checked constructor cannot reject it.
            ErrorRate::try_new(10f64.powf(-phred / 10.0))
                .expect("a Phred between 10 and 50 is a probability in [0, 1]")
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The ladder's shape is a correctness parameter of the scan, not a
    /// convenience: 161 rungs at quarter-Phred spacing is what
    /// `spec/parameter_prepass.md` §3 argues is finer than a caller can feel, and
    /// the whole "no refinement stage" decision rests on it.
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

    /// Ascending in Phred means **descending** in probability, and the ratio
    /// between neighbours is constant — that is what makes "one rung" a meaningful
    /// unit of distance for the fits that report movement in rungs.
    #[test]
    fn the_error_rate_ladder_rungs_are_a_constant_ratio_apart() {
        let ladder = error_rate_ladder();
        let expected_ratio = 10f64.powf(0.025);

        for pair in ladder.windows(2) {
            let (coarser, finer) = (pair[0].get(), pair[1].get());
            assert!(
                finer < coarser,
                "rates descend as Phred ascends: {finer} should be below {coarser}"
            );
            let ratio = coarser / finer;
            assert!(
                (ratio - expected_ratio).abs() < 1e-12,
                "adjacent rungs differ by 10^0.025, got {ratio}"
            );
        }
    }

    /// A window is 100 kb. Stated as a test because the number is load-bearing:
    /// the runs model's noise floor is a function of how many windows a genome has
    /// (`spec/parameter_prepass_generic.md` §6.1), so changing this changes the
    /// resolution of every reported inbreeding coefficient.
    #[test]
    fn the_inbreeding_window_is_a_hundred_kilobases() {
        assert_eq!(INBREEDING_WINDOW_BP.get(), 100_000);
    }
}
