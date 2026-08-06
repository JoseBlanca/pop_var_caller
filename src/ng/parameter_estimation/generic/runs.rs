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

/// The starting points the runs model climbs from.
///
/// **They must disagree about how far apart the two states are, not only about how much
/// of the genome is inside a run — and that is the whole content of this type.** Starts
/// that differ only in the second are not a spread: they all miss a genome whose real
/// states sit close together, in the same way, and "keep the best-scoring" then has
/// nothing better to pick from. Measured: nine starts spanning the separation return
/// `F` = 0.2634 where five sharing one separation guess return `F` = 0.0000 —
/// converged, silent, and on the same genome, 29% of which is covered by runs.
///
/// The defaults are three separations crossed with three inside fractions, and they
/// cost seconds on 8,000 windows.
#[derive(Clone, PartialEq, Debug)]
pub struct RunsModelStarts {
    /// The inside state's heterozygote rate as a fraction of the outside state's — how
    /// far apart the two states are guessed to be.
    pub separations: SmallVec<[f64; 3]>,
    /// The fraction of the genome each start begins by assuming lies inside a run.
    pub implied_f: SmallVec<[f64; 3]>,
}

impl Default for RunsModelStarts {
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
    /// The heterozygote rate **inside** a run. Fitted rather than fixed at zero: at
    /// zero, one collapsed paralog inside a run costs the whole
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
    /// rejected. Sorted by score, best first.
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
}
