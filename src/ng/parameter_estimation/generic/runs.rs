//! The inbreeding coefficient: a two-state hidden Markov model over windows, deciding
//! which of them lie in a **run of homozygosity** — a stretch of genome where both
//! copies descend from one recent ancestor, and which therefore carries almost no
//! heterozygotes. The rest of the genome carries them at the sample's own rate.
//!
//! The chain walks the windows of a contig deciding which state each is in, and the
//! coefficient is the share of the analysable genome the inside state claims —
//! weighted by reference positions rather than by loci, so a window dense in widened
//! indel loci is not under-weighted (`spec/parameter_prepass_generic.md` §6).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §5.3. The fit itself lands
//! in Milestone E; the types below are Milestone A.

use smallvec::SmallVec;

use crate::ng::parameter_estimation::fitting::FitTermination;

/// Fewest windows the runs model will accept.
///
/// **Not the same kind of floor as
/// [`MIN_SITES_TO_FIT`](super::MIN_SITES_TO_FIT).** That one guards against a rate
/// fitted from too little; this one guards against an answer that is entirely the
/// estimator's own noise. A genome generated with **no runs at all** returned `F`
/// averaging 0.23 at 1,200 windows, and 0.84 on one seed of eight
/// (`spec/parameter_prepass_generic.md` §6.1). A tomato genome is 8,004 windows and a
/// human 31,000, so no real run comes near this — but a development fixture or a
/// region-restricted run does, and the number it would produce looks like any other.
/// Fail rather than emit.
pub const MIN_WINDOWS_TO_FIT_INBREEDING: usize = 3_000;

/// The starting points the runs model climbs from.
///
/// **They must disagree about how far apart the two states are, not only about how much
/// of the genome is inside a run — and that is the whole content of this type.** Starts
/// that differ only in the second are not a spread: they all miss a genome whose real
/// states sit close together, in the same way, and "keep the best-scoring" then has
/// nothing better to pick from. Measured on a genome 26% covered by runs: nine starts
/// spanning the separation return `F` = 0.2634 against a realised 0.2629, where five
/// sharing one separation guess return `F` = 0.0000 — converged, and silent
/// (`research/parameter_estimator_experiments_2026-08-06.md` §3.4).
///
/// The default is three separations crossed with three inside fractions, and it costs
/// seconds on 8,000 windows.
///
/// **A caller who replaces the default must keep the separations distinct.** Nothing
/// here enforces it — the fit is what rejects a start set that leaves one state empty,
/// with `InbreedingStatesNotSeparated` rather than a zero.
#[derive(Clone, PartialEq, Debug)]
pub struct RunsModelStarts {
    /// The inside state's heterozygote rate **as a fraction of the outside state's**.
    ///
    /// **A smaller number means the two states are guessed further apart**, which is the
    /// opposite of what the word "separation" suggests: 0.05 guesses the inside state at
    /// a twentieth of the outside rate, and 0.75 at three quarters of it. Getting this
    /// backwards is how a start set ends up spanning nothing.
    ///
    /// Defaults to `[0.05, 1/3, 0.75]`.
    pub separations: SmallVec<[f64; 3]>,
    /// The fraction of the genome each start begins by assuming lies inside a run.
    /// Defaults to `[0.05, 0.5, 0.75]`.
    pub implied_f: SmallVec<[f64; 3]>,
}

impl Default for RunsModelStarts {
    /// Three separations crossed with three inside fractions — nine starts. The values
    /// are on the two fields; the property that matters is that the *separations*
    /// differ, not that there are nine.
    fn default() -> Self {
        Self {
            separations: SmallVec::from_slice(&[0.05, 1.0 / 3.0, 0.75]),
            implied_f: SmallVec::from_slice(&[0.05, 0.5, 0.75]),
        }
    }
}

/// One starting point's result, kept so that a search which found nothing can be told
/// apart from a genome that has nothing to find.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct StartOutcome {
    pub separation: f64,
    pub implied_f: f64,
    pub inbreeding: f64,
    pub log_likelihood: f64,
}

/// What the runs model fitted alongside `F`, and how it terminated.
///
/// Emitted because every one of these is a number someone will want when `F` looks
/// wrong.
#[derive(Clone, PartialEq, Debug)]
pub struct RunsModelFit {
    /// The heterozygote rate outside a run — the sample's ordinary heterozygosity.
    pub outside_het: f64,
    /// The homozygous-non-reference rate outside a run.
    pub outside_hom_alt: f64,
    /// The heterozygote rate **inside** a run — the small residue of apparent
    /// heterozygotes a truly autozygous stretch still shows, from collapsed paralogs and
    /// mismapping.
    ///
    /// **A fitted rate, and "floor" names its role rather than its arithmetic**: it is
    /// what stops the inside state charging an impossible cost. Fitted rather than fixed
    /// at zero, because at zero one collapsed paralog inside a run costs the whole
    /// heterozygote-against-homozygous-reference ratio at that site — about 125 for one
    /// alternative read of three, and past 10³⁵ for fifteen of thirty.
    pub inside_het_floor: f64,
    /// The homozygous-non-reference rate inside a run.
    pub inside_hom_alt: f64,
    /// The fitted per-base chance of entering a run.
    pub enter_run_per_base: f64,
    /// The fitted per-base chance of leaving one.
    pub leave_run_per_base: f64,
    pub termination: FitTermination,

    /// **How the search went, and the only thing separating a real `F` = 0 from a
    /// failed one.** Both leave the inside state empty and its frequencies at their
    /// starting values; only the scores say whether a better answer was looked for and
    /// rejected. Sorted by score, best first — read the winner through
    /// [`Self::best_start`] rather than indexing, so the ordering is stated in one place.
    pub starts_tried: SmallVec<[StartOutcome; 9]>,

    /// The noise floor at this run's window count — what `F` comes back as on a genome
    /// with no runs at all. About 0.01 at tomato's 8,004 windows and 0.003 at a human
    /// genome's 31,000. **An `F` below it is *nothing detected*, and a consumer that
    /// cannot see this number cannot tell that from a small autozygous fraction.**
    pub resolution: f64,

    /// Windows whose posterior landed between 0.01 and 0.99 — the ones the chain rather
    /// than their own reads decided. **Zero at 100 kb**, which is the measurement saying
    /// the transitions changed no window's classification. Non-zero is not a fault; it
    /// is the chain earning its keep.
    pub undecided_windows: u32,
}

impl RunsModelFit {
    /// The starting point that scored best — the one `F` was taken from.
    ///
    /// The single reader of `starts_tried`'s "best first" ordering, so that the ordering
    /// is a property of this file rather than of every call site that indexes. `None`
    /// only if no start was tried at all, which the fit does not produce.
    #[must_use]
    pub fn best_start(&self) -> Option<&StartOutcome> {
        self.starts_tried.first()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The default starts span the state separation, and that is the property the
    /// type exists for.** Three distinct separations is what makes them a spread; three
    /// distinct inside fractions sharing one separation would not be, and the
    /// difference is `F` = 0.2634 against a confident `F` = 0.0000 on the same genome.
    #[test]
    fn the_default_starts_disagree_about_the_state_separation_not_only_about_f() {
        let starts = RunsModelStarts::default();

        assert_eq!(starts.separations.len(), 3);
        assert_eq!(starts.implied_f.len(), 3);
        for pair in starts.separations.windows(2) {
            assert!(
                pair[0] < pair[1],
                "the separations must differ, and ascend: {:?}",
                starts.separations
            );
        }
        for pair in starts.implied_f.windows(2) {
            assert!(pair[0] < pair[1], "{:?}", starts.implied_f);
        }
    }

    /// Every start is a real guess: a separation is a fraction of the outside rate and
    /// an inside fraction is a fraction of the genome, so both live strictly inside
    /// `(0, 1)`. A separation of 1 would make the two states identical, and one of 0
    /// would fix the inside state's heterozygote rate at zero — the tie the fit is
    /// written to avoid.
    #[test]
    fn every_default_start_lies_strictly_inside_the_unit_interval() {
        let starts = RunsModelStarts::default();

        for guess in starts.separations.iter().chain(starts.implied_f.iter()) {
            assert!(*guess > 0.0 && *guess < 1.0, "{guess} is not a real guess");
        }
        assert_eq!(starts.separations.len() * starts.implied_f.len(), 9);
    }

    /// **The two fields hold the same type and the same shape, so nothing but their
    /// values tells them apart.** Both are three ascending fractions in `(0, 1)`, so
    /// every property asserted above survives swapping them — and a swap puts the spread
    /// on the wrong axis, which is precisely the `F` = 0.0000 this type exists to
    /// prevent. The values are therefore pinned.
    #[test]
    fn the_default_starts_are_the_values_the_design_specifies_on_each_axis() {
        let starts = RunsModelStarts::default();

        assert_eq!(starts.separations.as_slice(), &[0.05, 1.0 / 3.0, 0.75]);
        assert_eq!(starts.implied_f.as_slice(), &[0.05, 0.5, 0.75]);
    }

    /// `starts_tried` is documented best-first, and `best_start` is the one reader of
    /// that ordering — so a call site cannot index the wrong end.
    #[test]
    fn best_start_reads_the_head_of_the_best_first_ordering() {
        let outcome = |inbreeding: f64, log_likelihood: f64| StartOutcome {
            separation: 1.0 / 3.0,
            implied_f: 0.5,
            inbreeding,
            log_likelihood,
        };
        let fit = RunsModelFit {
            outside_het: 0.0105,
            outside_hom_alt: 0.0010,
            inside_het_floor: 0.0008,
            inside_hom_alt: 0.0009,
            enter_run_per_base: 1e-7,
            leave_run_per_base: 3e-7,
            termination: FitTermination {
                iterations: 12,
                converged: true,
            },
            starts_tried: SmallVec::from_slice(&[
                outcome(0.2634, -1.20e9),
                outcome(0.0000, -1.41e9),
            ]),
            resolution: 0.01,
            undecided_windows: 0,
        };

        assert_eq!(fit.best_start().map(|s| s.inbreeding), Some(0.2634));
        assert!(
            RunsModelFit {
                starts_tried: SmallVec::new(),
                ..fit
            }
            .best_start()
            .is_none()
        );
    }

    /// The runs model's floor is stated where the model lives, not beside the
    /// site-count floor it is unlike: that one guards against a rate fitted from too
    /// little, this one against an answer that is entirely the estimator's own noise.
    #[test]
    fn the_window_floor_is_the_measured_one() {
        assert_eq!(MIN_WINDOWS_TO_FIT_INBREEDING, 3_000);
    }
}
