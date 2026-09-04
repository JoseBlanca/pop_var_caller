//! Re-fitting one repeat tract's slippage numbers from its own reads — the arithmetic of the
//! calling loop's outer round (`doc/devel/ng/spec/calling_em_loop.md` §5.1).
//!
//! **This module holds the counts and the formulas; the round itself lives in the loop's
//! driver** (`summarise_condition`), because attributing a read needs the genotype posteriors
//! and the cached emissions, and both are the driver's to hand over. What lives here is what
//! can be stated and tested without a locus: how the pooled counts become three re-fitted
//! numbers, how hard they are pulled back toward the stratum's frozen values, and when the
//! rounds have stopped moving.
//!
//! # The granularity: one pooled set of counts per locus, applied to every cell
//!
//! **Production re-fits one set of numbers per locus, pooled over every read group** — its
//! attribution accumulates every sample's reads into a single `LocusSlipFit`
//! (`src/ssr/cohort/em.rs:1200`, one struct built before the sample loop and never keyed by
//! group), and each round re-fits one shape and one level multiplier from it
//! (`em.rs:574-575`). The per-group structure lives entirely in the *frozen* inputs: the
//! group's own level line is scaled by the one locus-wide multiplier. ng mirrors that shape.
//! The counts below are pooled over the locus's read groups; the level moves as one multiplier
//! on every cell's frozen level; and where production pulls its one shape toward a per-period
//! parent (`theta0`), ng pulls **each `(read group, candidate)` cell's shape toward that
//! cell's own frozen stratum values** — ng has no per-period parent shape, and the stratum's
//! frozen values are what spec §5.1 names as the pull-back target. At zero slips every cell
//! therefore collapses to its own frozen numbers, exactly as production collapses to its
//! parent shape (`src/ssr/cohort/stutter.rs:188`).
//!
//! # What differs from production, named
//!
//! - **Attribution is posterior-weighted, not hard.** Production attributes each read to the
//!   *called* genotype's nearest allele (`em.rs:1192`); ng splits each read's weight across
//!   the genotype posteriors, and within a genotype by the responsibility
//!   `copy share × cached emission` — the spec's explicit ruling (§5.1), and the refinement
//!   production's own comment defers (`em.rs:18`). The counts here are therefore fractional
//!   `f64` weights where production's are whole reads.
//! - **A part-repeat difference is out of the re-fit.** Production floors a read's length
//!   into whole units (`obs.len() / period`) and fits it anyway; ng excludes any
//!   `(read, attributed allele)` pair whose length difference is not a whole number of motif
//!   units, because the part-repeat shares are placeholders the re-fit must not fit
//!   (spec §5.1 — *"hold the count at three"*).

use crate::ng::calling::inference::SlippageRefitConfig;
use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;

/// The largest slip, in whole motif units, the fall-off histogram counts — production's
/// `MAX_SLIP` (`src/ssr/cohort/param_estimation.rs:21`).
///
/// **A slip past it still counts as slipped and leaves the shape alone**, mirroring
/// production's asymmetry (`param_estimation.rs:222-231` drops it from the histogram;
/// `em.rs:1226` counts it toward the level regardless): the read clearly slipped, but a
/// magnitude the geometric tail is truncated at would drag the mean it is fitted from.
pub(crate) const MAX_SLIP_UNITS: u64 = 10;

/// The most the locus's own reads may scale a cell's frozen slippage level — production's
/// `LEVEL_MULT_MAX` (`src/ssr/cohort/em.rs:1114`).
pub(crate) const LEVEL_MULTIPLIER_CEILING: f64 = 10.0;

/// Where the raw fall-off estimate is clamped before the pull-back blends it — production's
/// bound (`src/ssr/cohort/stutter.rs:218`). The blend below is convex, so the blended value
/// stays inside the range whenever the frozen value does.
const FALL_OFF_ESTIMATE_CEILING: f64 = 0.999;

/// **The locus's pooled slip evidence**: every read's posterior-weighted attribution, summed
/// over the whole locus in the driver's fixed row-and-observation order.
///
/// The weights are fractional — a read split across two genotypes contributes to both — so
/// every field is an `f64` accumulated in one fixed order, which is what keeps a re-fit
/// identical at any worker count (spec §8).
#[derive(Debug, Default, Clone, PartialEq)]
pub(crate) struct PooledSlipCounts {
    /// Total weight attributed to any allele — every read the re-fit saw, slipped or not.
    pub(crate) attributed: f64,
    /// Weight attributed at a length other than its allele's — the level's numerator.
    pub(crate) slipped: f64,
    /// Weight × the attributed cell's **frozen** level, summed — how many slips the frozen
    /// numbers predicted for these reads, the level's denominator. Production's
    /// `expected_slipped` (`em.rs:1232`), which reads the raw group line with no multiplier —
    /// so each round's multiplier replaces the last rather than compounding onto it.
    pub(crate) expected_slipped: f64,
    /// Slipped weight within [`MAX_SLIP_UNITS`] — the shape's `n`.
    pub(crate) in_histogram: f64,
    /// Of that, the weight that slipped **shorter** — the direction split's numerator.
    pub(crate) shorter: f64,
    /// Σ weight × |slip| in whole units, over the histogram — what the fall-off's mean
    /// magnitude is read from.
    pub(crate) magnitude: f64,
}

impl PooledSlipCounts {
    /// Add one read's weight, attributed to an allele `slip_units` whole motif units away,
    /// whose cell's frozen slippage level is `frozen_level`.
    ///
    /// The caller has already excluded what does not belong here: a part-repeat difference,
    /// a read that ran out inside the tract, and a cell the fit has no numbers for.
    pub(crate) fn add(&mut self, slip_units: i64, weight: f64, frozen_level: f64) {
        self.attributed += weight;
        self.expected_slipped += weight * frozen_level;
        if slip_units == 0 {
            return;
        }
        self.slipped += weight;
        let size = slip_units.unsigned_abs();
        if size > MAX_SLIP_UNITS {
            // Production's asymmetry, kept deliberately: an over-cap slip counts toward the
            // level and not toward the shape (`param_estimation.rs:222-231`).
            return;
        }
        self.in_histogram += weight;
        if slip_units < 0 {
            self.shorter += weight;
        }
        self.magnitude += weight * size as f64;
    }
}

/// **How far the locus's reads scale every cell's frozen level** — production's formula
/// verbatim, with the pull-back a setting instead of a constant:
/// `clamp((slipped + P) / (expected_slipped + P), 0, 10)`
/// (`src/ssr/cohort/em.rs:1329-1331`, `P` = 20 slipped reads there).
///
/// The pull-back's target is a multiplier of **one** — the frozen level unchanged — and the
/// frozen level at each cell is already the fitted curve's value read off at that cell
/// (`str_slippage_level_curve.md`; the gather stores what `StratumFits::at` blended), so
/// pulling the multiplier to one is pulling the level to the curve, which is what
/// [`DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS`](super::DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS)
/// documents.
///
/// **A locus with no evidence and no pull-back leaves the level alone**: at a denominator of
/// zero the ratio is no estimate at all, and one is the multiplier that changes nothing.
pub(crate) fn level_multiplier(counts: &PooledSlipCounts, pull_back_slipped_reads: f64) -> f64 {
    let denominator = counts.expected_slipped + pull_back_slipped_reads;
    if denominator <= 0.0 {
        return 1.0;
    }
    ((counts.slipped + pull_back_slipped_reads) / denominator).clamp(0.0, LEVEL_MULTIPLIER_CEILING)
}

/// One cell's three re-fitted numbers, from the locus's pooled counts and that cell's frozen
/// values.
///
/// The formulas are production's (`src/ssr/cohort/stutter.rs:184-223`), with the pseudo-count
/// a setting and the pull-back target this cell's frozen stratum values (see the module note
/// for why the target differs from production's per-period parent):
///
/// - **level** — the frozen level scaled by the one locus-wide multiplier, clamped into
///   `[0, 1]` as production clamps its per-read level (`em.rs:385-388`);
/// - **direction split** — `(shorter + S·frozen) / (n + S)`, the Beta posterior mean with `S`
///   pseudo-slips at the frozen split;
/// - **fall-off** — the geometric continuation probability whose mean magnitude matches the
///   histogram's, `clamp((mean − 1) / mean, 0, 0.999)`, blended `(n·estimate + S·frozen) / (n + S)`.
///
/// **A locus whose histogram is empty keeps the frozen shape verbatim**, as production returns
/// its prior at `total == 0` (`stutter.rs:188-190`) — which also answers the free setting's
/// `0/0`: no slips and no pseudo-counts is no estimate, not a zero.
pub(crate) fn refit_cell(
    frozen: &Slippage,
    counts: &PooledSlipCounts,
    multiplier: f64,
    shape_pull_back_pseudocounts: f64,
) -> Slippage {
    let level = (frozen.level * multiplier).clamp(0.0, 1.0);
    let n = counts.in_histogram;
    let (shorter_share, fall_off) = if n > 0.0 {
        let pull = shape_pull_back_pseudocounts;
        let shorter_share = (counts.shorter + pull * frozen.shorter_share) / (n + pull);
        let mean_magnitude = counts.magnitude / n;
        let estimate =
            ((mean_magnitude - 1.0) / mean_magnitude).clamp(0.0, FALL_OFF_ESTIMATE_CEILING);
        let fall_off = (n * estimate + pull * frozen.fall_off) / (n + pull);
        (shorter_share, fall_off)
    } else {
        (frozen.shorter_share, frozen.fall_off)
    };
    Slippage {
        level,
        shorter_share,
        fall_off,
    }
}

/// Every cell's re-fitted numbers, into a buffer the driver reuses across rounds.
///
/// **A cell the fit has no frozen numbers for stays `None`** — it keeps the shipped constants
/// it was scored under, because there is nothing of its own to pull a re-fit toward, and its
/// reads were left out of the counts for the same reason.
pub(crate) fn refit_cells(
    frozen_of_each_cell: &[Option<Slippage>],
    counts: &PooledSlipCounts,
    config: &SlippageRefitConfig,
    refitted: &mut Vec<Option<Slippage>>,
) {
    let multiplier = level_multiplier(counts, config.level_pull_back_slipped_reads);
    refitted.clear();
    refitted.extend(frozen_of_each_cell.iter().map(|cell| {
        cell.as_ref().map(|frozen| {
            refit_cell(
                frozen,
                counts,
                multiplier,
                config.direction_and_fall_off_pull_back_pseudocounts,
            )
        })
    }));
}

/// The largest absolute move any re-fitted number made — what the round's stopping rule
/// compares against
/// [`round_convergence_threshold`](super::SlippageRefitConfig::round_convergence_threshold).
///
/// Absolute movement on every number, as production compares its four
/// (`src/ssr/cohort/em.rs:576`, `shapes_close` at `em.rs:1334`); ng's re-fitted numbers live
/// per cell, so *"every re-fitted number moves less than a threshold"* (spec §6) is read over
/// every cell's three.
///
/// # Panics
///
/// If the two slices disagree about the cells — length, or which cells carry numbers. Both
/// come from the same gather over the same locus, so a disagreement is two loci mixed up.
pub(crate) fn largest_movement(current: &[Option<Slippage>], refitted: &[Option<Slippage>]) -> f64 {
    assert_eq!(
        current.len(),
        refitted.len(),
        "the two rounds' cell tables cover {} and {} cells, so they describe different loci",
        current.len(),
        refitted.len()
    );
    let mut largest: f64 = 0.0;
    for (before, after) in current.iter().zip(refitted) {
        match (before, after) {
            (None, None) => {}
            (Some(before), Some(after)) => {
                largest = largest
                    .max((after.level - before.level).abs())
                    .max((after.shorter_share - before.shorter_share).abs())
                    .max((after.fall_off - before.fall_off).abs());
            }
            (None, Some(_)) | (Some(_), None) => panic!(
                "a cell carries fitted numbers in one round and none in the other; whether a \
                 cell has frozen numbers is a property of the gather, so the two rounds \
                 describe different loci"
            ),
        }
    }
    largest
}

/// **Every sample's cached emissions for one locus, per `(row, observation, candidate)`** —
/// copied out of the row scratch as the table build fills it, so the re-fit's attribution can
/// read the same numbers the genotype likelihoods were assembled from.
///
/// It exists because the worker's row scratch is one buffer reused across rows: after the
/// table build only the last sample's emissions survive there
/// (`CallingScratch::ssr_row`), and the spec's responsibility split reads *"that read's cached
/// emission for the allele"* (§5.1). The copy costs one `memcpy` per row per round and is only
/// made where a re-fit round asked for it — the frozen path never touches this type.
#[derive(Debug, Default)]
pub(crate) struct RefitEmissionCache {
    /// Every row's emissions, concatenated in row order; each row is
    /// `observations × candidates`, observation-major — the row scratch's own layout.
    emissions: Vec<f64>,
    /// Where each row starts in [`Self::emissions`].
    row_starts: Vec<usize>,
    /// The stride: how many candidates each observation's block holds.
    candidates: usize,
}

impl RefitEmissionCache {
    /// Clear the cache for one locus's build.
    pub(crate) fn begin_locus(&mut self, candidates: usize) {
        self.emissions.clear();
        self.row_starts.clear();
        self.candidates = candidates;
    }

    /// Append one row's freshly filled emissions, in the build's own row order.
    pub(crate) fn push_row(&mut self, row_emissions: &[f64]) {
        self.row_starts.push(self.emissions.len());
        self.emissions.extend_from_slice(row_emissions);
    }

    /// One cached emission — the observation axis is every observation of the sample,
    /// partials included, exactly as the row scratch indexes it.
    ///
    /// # Panics
    ///
    /// On a row, observation or candidate past what the build pushed — held in release, as
    /// every scratch index on this path is (spec §8).
    pub(crate) fn emission_at(&self, row: usize, observation: usize, candidate: usize) -> f64 {
        assert!(
            candidate < self.candidates,
            "candidate {candidate} is past the {} this cache was built over",
            self.candidates
        );
        let start = self.row_starts[row];
        let end = self
            .row_starts
            .get(row + 1)
            .copied()
            .unwrap_or(self.emissions.len());
        let slot = start + observation * self.candidates + candidate;
        assert!(
            slot < end,
            "observation {observation} of row {row} addresses slot {slot}, past that row's end \
             at {end} — the cache and the evidence describe different loci"
        );
        self.emissions[slot]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn a_frozen_cell() -> Slippage {
        Slippage {
            level: 0.04,
            shorter_share: 0.83,
            fall_off: 0.25,
        }
    }

    fn shipped_config_with_rounds() -> SlippageRefitConfig {
        SlippageRefitConfig {
            max_rounds: 3,
            ..SlippageRefitConfig::DEFAULT
        }
    }

    /// **The level formula is production's with the pull-back a setting**: 10 slipped reads
    /// where the frozen numbers expected 1.6 gives `(10 + 20) / (1.6 + 20) = 1.3888…` at the
    /// shipped pull-back, so a frozen level of 0.04 becomes `0.0555…` — moved toward the
    /// reads' own rate of 0.25, and held far short of it by the 20 pseudo-slips.
    #[test]
    fn the_level_moves_toward_the_reads_and_is_held_by_the_pull_back() {
        let mut counts = PooledSlipCounts::default();
        // Forty reads at a cell whose frozen level is 0.04: ten slipped one unit shorter.
        for _ in 0..30 {
            counts.add(0, 1.0, 0.04);
        }
        for _ in 0..10 {
            counts.add(-1, 1.0, 0.04);
        }
        let multiplier = level_multiplier(&counts, 20.0);
        assert!((multiplier - 30.0 / 21.6).abs() < 1e-12, "{multiplier}");
        let refitted = refit_cell(&a_frozen_cell(), &counts, multiplier, 50.0);
        assert!(
            (refitted.level - 0.04 * (30.0 / 21.6)).abs() < 1e-12,
            "{}",
            refitted.level
        );
        // And the shape: ten shorter slips against 50 pseudo-counts at 0.83 —
        // (10 + 41.5) / 60; every slip was one unit, so the raw fall-off estimate is 0 and
        // the blend is (0 + 50 × 0.25) / 60.
        assert!((refitted.shorter_share - 51.5 / 60.0).abs() < 1e-12);
        assert!((refitted.fall_off - 12.5 / 60.0).abs() < 1e-12);
    }

    /// **Zero pull-back is the free setting** — the numbers go to the reads' own rate
    /// (HipSTR's behaviour, spec §5.1): the multiplier is `slipped / expected`, the direction
    /// split the observed share, the fall-off the raw geometric estimate.
    #[test]
    fn zero_pull_back_hands_the_numbers_to_the_reads() {
        let mut counts = PooledSlipCounts::default();
        for _ in 0..30 {
            counts.add(0, 1.0, 0.04);
        }
        // Ten slips: eight one unit shorter, two two units longer — mean magnitude 1.2.
        for _ in 0..8 {
            counts.add(-1, 1.0, 0.04);
        }
        counts.add(2, 1.0, 0.04);
        counts.add(2, 1.0, 0.04);
        let multiplier = level_multiplier(&counts, 0.0);
        assert!((multiplier - 10.0 / 1.6).abs() < 1e-12);
        let refitted = refit_cell(&a_frozen_cell(), &counts, multiplier, 0.0);
        // 0.04 × 6.25 = 0.25: the reads' own slip rate, 10 slipped of 40.
        assert!((refitted.level - 0.25).abs() < 1e-12);
        assert!((refitted.shorter_share - 0.8).abs() < 1e-12);
        // mean magnitude 1.2 → (1.2 − 1) / 1.2.
        assert!((refitted.fall_off - 0.2 / 1.2).abs() < 1e-12);
    }

    /// **No slips leaves the shape frozen verbatim** (production returns its prior at
    /// `total == 0`), and pulls the level only as far as the evidence of clean reads warrants:
    /// ten clean reads at a frozen level of 0.04 expected 0.4 slips, so the multiplier is
    /// `20 / 20.4` and the level eases to `0.0392…` rather than snapping to anything.
    #[test]
    fn a_locus_with_no_slips_keeps_the_frozen_shape() {
        let mut counts = PooledSlipCounts::default();
        for _ in 0..10 {
            counts.add(0, 1.0, 0.04);
        }
        let multiplier = level_multiplier(&counts, 20.0);
        assert!((multiplier - 20.0 / 20.4).abs() < 1e-12);
        let refitted = refit_cell(&a_frozen_cell(), &counts, multiplier, 50.0);
        assert_eq!(refitted.shorter_share, 0.83);
        assert_eq!(refitted.fall_off, 0.25);
        assert!((refitted.level - 0.04 * 20.0 / 20.4).abs() < 1e-12);
    }

    /// **No evidence at all changes nothing at any pull-back** — the multiplier is one, the
    /// shape is frozen, and the free setting's `0/0` is *no estimate* rather than a zero.
    #[test]
    fn no_evidence_at_all_is_no_estimate() {
        let counts = PooledSlipCounts::default();
        for pull in [20.0, 0.0] {
            let multiplier = level_multiplier(&counts, pull);
            assert_eq!(multiplier, 1.0, "at a pull-back of {pull}");
        }
        let refitted = refit_cell(&a_frozen_cell(), &counts, 1.0, 0.0);
        assert_eq!(refitted, a_frozen_cell());
    }

    /// **The multiplier is capped at production's ceiling of 10**, so a locus of nothing but
    /// slips cannot scale a cell's level without bound — and the level itself is clamped into
    /// `[0, 1]` even where the frozen value times ten would leave it.
    #[test]
    fn the_multiplier_and_the_level_are_both_capped() {
        let mut counts = PooledSlipCounts::default();
        for _ in 0..100 {
            counts.add(-1, 1.0, 0.0005);
        }
        assert_eq!(level_multiplier(&counts, 0.0), LEVEL_MULTIPLIER_CEILING);
        let heavy = Slippage {
            level: 0.2,
            ..a_frozen_cell()
        };
        let refitted = refit_cell(&heavy, &counts, LEVEL_MULTIPLIER_CEILING, 0.0);
        assert_eq!(refitted.level, 1.0);
    }

    /// **A slip past ten whole units counts as slipped and stays out of the shape** —
    /// production's asymmetry, mirrored on purpose.
    #[test]
    fn an_over_cap_slip_reaches_the_level_and_not_the_shape() {
        let mut counts = PooledSlipCounts::default();
        counts.add(-11, 1.0, 0.04);
        assert_eq!(counts.slipped, 1.0);
        assert_eq!(counts.in_histogram, 0.0);
        assert_eq!(counts.shorter, 0.0);
        assert_eq!(counts.magnitude, 0.0);
    }

    /// A cell with no frozen numbers stays `None` through a re-fit, and the movement test
    /// reads it as no movement.
    #[test]
    fn a_cell_the_fit_never_reached_is_left_alone() {
        let mut counts = PooledSlipCounts::default();
        counts.add(-1, 4.0, 0.04);
        let frozen = [Some(a_frozen_cell()), None];
        let mut refitted = Vec::new();
        refit_cells(
            &frozen,
            &counts,
            &shipped_config_with_rounds(),
            &mut refitted,
        );
        assert_eq!(refitted.len(), 2);
        assert!(refitted[0].is_some());
        assert!(refitted[1].is_none());
        assert!(largest_movement(&frozen, &refitted) > 0.0);
        assert_eq!(largest_movement(&refitted, &refitted), 0.0);
    }

    /// The movement is the largest absolute change over every cell's three numbers —
    /// production's rule (`shapes_close`, absolute and simultaneous), read over the cells.
    #[test]
    fn the_movement_is_the_largest_absolute_change_of_any_number() {
        let before = [Some(a_frozen_cell())];
        let after = [Some(Slippage {
            level: 0.041,
            shorter_share: 0.85,
            fall_off: 0.24,
        })];
        let movement = largest_movement(&before, &after);
        assert!((movement - 0.02).abs() < 1e-12, "{movement}");
    }

    /// The emission cache hands back what each row's build pushed, at the row scratch's own
    /// `(observation, candidate)` stride.
    #[test]
    fn the_emission_cache_reads_back_what_the_build_pushed() {
        let mut cache = RefitEmissionCache::default();
        cache.begin_locus(2);
        cache.push_row(&[0.1, 0.2, 0.3, 0.4]); // two observations × two candidates
        cache.push_row(&[0.5, 0.6]); // one observation
        assert_eq!(cache.emission_at(0, 0, 1), 0.2);
        assert_eq!(cache.emission_at(0, 1, 0), 0.3);
        assert_eq!(cache.emission_at(1, 0, 1), 0.6);
    }

    /// An observation past its own row's end is refused rather than read from the next row's
    /// block — held in release, like every scratch index on this path.
    #[test]
    #[should_panic(expected = "past that row's end")]
    fn an_observation_past_its_rows_end_is_refused() {
        let mut cache = RefitEmissionCache::default();
        cache.begin_locus(2);
        cache.push_row(&[0.1, 0.2]);
        cache.push_row(&[0.3, 0.4]);
        let _ = cache.emission_at(0, 1, 0);
    }
}
