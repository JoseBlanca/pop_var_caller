//! **Each read group's stutter multiplier for one motif length, with the classes' rates held**
//! (spec `str_library_stutter.md` §3.4–§3.6).
//!
//! A read group's reads slip at its class's rate times the read group's multiplier, so a read
//! group whose polymerase or PCR slips twice as often as typical has a multiplier of two. With the
//! class rates `Lₛ` held, read group `r`'s log-likelihood in its multiplier `m` is
//!
//! `Σₛ [ Sᵣₛ · ln(m Lₛ) + Nᵣₛ · ln(1 − m Lₛ) ] + P · (ln m − m)`,
//!
//! `Sᵣₛ` and `Nᵣₛ` the read group's expected slipped and exactly-read reads in class `s` from the
//! expectation step, and the last term the pull toward one: `P` pseudo-reads that slipped as often
//! as a typical read group's. It is concave in `ln m`, so it is climbed by Newton's method on
//! `ln m`.
//!
//! **What fixes the scale.** `m · Lₛ` is unchanged when every multiplier is doubled and every class
//! rate halved, so a motif length's multipliers are scaled so that their geometric mean, each read
//! group weighted by its reads, is one ([`scale_to_a_typical_read_group`]).
//!
//! Nothing here reads the evidence: these are the functions the penalised fit calls between drawing
//! its curves (`penalised_fit.rs`).

use super::expected_counts::ReadGroupSlips;
use super::*;

/// **The pull toward a multiplier of one, in pseudo-reads** (spec §3.5): the per-locus slip refit's
/// own pull (`DEFAULT_LEVEL_PULL_BACK_SLIPPED_READS`). Soft; to be set from the stuttering
/// simulation at 3 reads a sample.
pub const PULL_TOWARD_ONE_SLIPPED_READS: f64 = 20.0;

/// **One class's evidence about one read group's multiplier**: the class's rate, held, and the read
/// group's expected exactly-read and slipped reads in it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClassEvidence {
    pub class_rate: f64,
    pub reads: ReadGroupSlips,
}

/// **The objective of the module documentation in `u = ln m`**, with its slope and its curvature
/// (the second derivative, negative): `(value, slope, curvature)`.
///
/// A class whose rate is zero carries no evidence about the multiplier and is skipped. At a
/// multiplier that lifts a class to `MULTIPLIED_LEVEL_CEILING` the objective is `−∞`; the solve
/// never steps there ([`largest_log_multiplier`]).
pub fn log_multiplier_objective(
    log_multiplier: f64,
    classes: &[ClassEvidence],
    pull: f64,
) -> (f64, f64, f64) {
    let multiplier = float::exp(log_multiplier);
    let (mut value, mut slope, mut curvature) = (
        pull * (log_multiplier - multiplier),
        pull * (1.0 - multiplier),
        -pull * multiplier,
    );
    for class in classes {
        if !(class.class_rate > 0.0) {
            continue;
        }
        let lifted = multiplier * class.class_rate;
        if lifted >= 1.0 {
            return (f64::NEG_INFINITY, f64::NAN, f64::NAN);
        }
        let ReadGroupSlips { exact, slipped } = class.reads;
        value += slipped * float::ln(lifted) + exact * float::ln_1p(-lifted);
        // d/du ln(1 − e^u L) = −e^u L / (1 − e^u L), and its derivative −e^u L / (1 − e^u L)².
        let odds = lifted / (1.0 - lifted);
        slope += slipped - exact * odds;
        curvature -= exact * odds / (1.0 - lifted);
    }
    (value, slope, curvature)
}

/// **The largest `ln m` the solve may reach**: the one that lifts the class of highest rate to
/// [`MULTIPLIED_LEVEL_CEILING`]. Infinite when no class has a rate.
pub fn largest_log_multiplier(classes: &[ClassEvidence]) -> f64 {
    let highest = classes
        .iter()
        .map(|class| class.class_rate)
        .filter(|rate| *rate > 0.0)
        .fold(0.0, f64::max);
    if highest > 0.0 {
        float::ln(MULTIPLIED_LEVEL_CEILING / highest)
    } else {
        f64::INFINITY
    }
}

/// **One read group's multiplier with the classes held**: the maximum of
/// [`log_multiplier_objective`], by Newton's method on `ln m` from `start`, each step halved until
/// the objective does not fall, and held below [`largest_log_multiplier`].
///
/// With no class evidence the pull alone holds it at one; with a pull of zero and no slipped read it
/// would run toward zero, so a run always pulls.
pub fn fit_multiplier(classes: &[ClassEvidence], pull: f64, start: f64) -> f64 {
    let ceiling = largest_log_multiplier(classes);
    // Just inside the ceiling, where every class's lifted rate is below one and the objective
    // finite.
    let below_ceiling = |u: f64| u.min(ceiling - 1e-9);
    let mut u = below_ceiling(float::ln(start));
    let (mut value, mut slope, mut curvature) = log_multiplier_objective(u, classes, pull);
    for _ in 0..100 {
        if !(curvature < 0.0) || slope.abs() < 1e-12 {
            break;
        }
        let mut step = -slope / curvature;
        let mut moved = false;
        for _ in 0..60 {
            let next = below_ceiling(u + step);
            let (next_value, next_slope, next_curvature) =
                log_multiplier_objective(next, classes, pull);
            if next_value >= value {
                moved = next != u;
                u = next;
                (value, slope, curvature) = (next_value, next_slope, next_curvature);
                break;
            }
            step *= 0.5;
        }
        if !moved {
            break;
        }
    }
    float::exp(u)
}

/// **Scales a motif length's multipliers so their geometric mean, each weighted by `weights`, is
/// one**, and returns the factor every multiplier was divided by (spec §2.2). A class rate times
/// that factor keeps every read group's slip rate where it was. Read groups with no weight do not
/// move the mean, and are scaled with the rest.
pub fn scale_to_a_typical_read_group(multipliers: &mut [f64], weights: &[f64]) -> f64 {
    assert_eq!(multipliers.len(), weights.len(), "one weight a read group");
    let total: f64 = weights.iter().sum();
    if !(total > 0.0) {
        return 1.0;
    }
    let mean_log = multipliers
        .iter()
        .zip(weights)
        .map(|(multiplier, weight)| weight * float::ln(*multiplier))
        .sum::<f64>()
        / total;
    let factor = float::exp(mean_log);
    for multiplier in multipliers.iter_mut() {
        *multiplier /= factor;
    }
    factor
}

/// **A multiplier's standard error on the log scale** (spec §3.6): `√((B + P·m) / H²)`, `H` minus the
/// objective's curvature at the fitted multiplier, pull included, `B` the sum over tracts of each
/// tract's contribution to the slope, squared, because the samples at one tract are not
/// independent, and `P·m` the pull's own curvature, which the stutter prior's errors also carry in
/// their middle term.
///
/// `by_tract` holds each tract's classes' evidence for this read group, one entry a tract with
/// reads from it. **`None` when the read group's slipped reads are fewer than the pull's
/// pseudo-reads**: its value then leans on the pull at least as much as on its reads.
pub fn log_multiplier_error(
    multiplier: f64,
    classes: &[ClassEvidence],
    by_tract: &[Vec<ClassEvidence>],
    pull: f64,
) -> Option<f64> {
    let slipped: f64 = classes.iter().map(|class| class.reads.slipped).sum();
    if slipped < pull {
        return None;
    }
    let u = float::ln(multiplier);
    let (_, _, curvature) = log_multiplier_objective(u, classes, pull);
    let information = -curvature;
    if !(information > 0.0) {
        return None;
    }
    let spread: f64 = by_tract
        .iter()
        .map(|tract| {
            let (_, slope, _) = log_multiplier_objective(u, tract, 0.0);
            slope * slope
        })
        .sum();
    Some(((spread + pull * multiplier) / (information * information)).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A small, seeded generator for drawn reads (splitmix64).
    struct Draw(u64);

    impl Draw {
        fn uniform(&mut self) -> f64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^= z >> 31;
            (z >> 11) as f64 / (1u64 << 53) as f64
        }

        fn normal(&mut self) -> f64 {
            let (a, b) = (self.uniform().max(1e-300), self.uniform());
            (-2.0 * a.ln()).sqrt() * (std::f64::consts::TAU * b).cos()
        }

        fn binomial(&mut self, trials: u32, chance: f64) -> u32 {
            (0..trials).filter(|_| self.uniform() < chance).count() as u32
        }
    }

    fn evidence(class_rate: f64, exact: f64, slipped: f64) -> ClassEvidence {
        ClassEvidence {
            class_rate,
            reads: ReadGroupSlips { exact, slipped },
        }
    }

    /// The class rates the drawn read groups are held at: from half a slip in a hundred reads to
    /// eight, as homopolymers run along the repeat count.
    fn class_rates() -> Vec<f64> {
        (0..20)
            .map(|class| 0.005 + 0.075 * class as f64 / 19.0)
            .collect()
    }

    /// One read group's reads drawn at `multiplier`: `tracts` tracts a class of `depth` reads each,
    /// every read slipping at its class's rate times the multiplier. Per tract, its classes'
    /// evidence (one class a tract); and over all of them, each class's.
    fn drawn_read_group(
        multiplier: f64,
        tracts: usize,
        depth: u32,
        draw: &mut Draw,
    ) -> (Vec<ClassEvidence>, Vec<Vec<ClassEvidence>>) {
        drawn_read_group_with_tracts_apart(multiplier, tracts, depth, 0.0, draw)
    }

    /// [`drawn_read_group`], each tract's rate its class's times a log-normal factor of mean one
    /// and log standard deviation `tract_spread`: tracts slipping more or less than their class,
    /// which is what makes the samples at one tract not independent.
    fn drawn_read_group_with_tracts_apart(
        multiplier: f64,
        tracts: usize,
        depth: u32,
        tract_spread: f64,
        draw: &mut Draw,
    ) -> (Vec<ClassEvidence>, Vec<Vec<ClassEvidence>>) {
        let mut classes = Vec::new();
        let mut by_tract = Vec::new();
        for rate in class_rates() {
            let (mut exact, mut slipped) = (0.0, 0.0);
            for _ in 0..tracts {
                let factor =
                    (tract_spread * draw.normal() - 0.5 * tract_spread * tract_spread).exp();
                let chance = (multiplier * rate * factor).min(0.95);
                let slips = f64::from(draw.binomial(depth, chance));
                let stays = f64::from(depth) - slips;
                exact += stays;
                slipped += slips;
                by_tract.push(vec![evidence(rate, stays, slips)]);
            }
            classes.push(evidence(rate, exact, slipped));
        }
        (classes, by_tract)
    }

    /// **The slope and the curvature are the objective's**, against central differences, at
    /// multipliers either side of one and with and without the pull.
    #[test]
    fn the_slope_and_curvature_are_the_objectives_own() {
        let classes = [
            evidence(0.01, 900.0, 12.0),
            evidence(0.05, 300.0, 30.0),
            evidence(0.2, 40.0, 15.0),
        ];
        for pull in [0.0, 20.0] {
            for u in [-0.7, 0.0, 0.4] {
                let step = 1e-5;
                let (_, slope, curvature) = log_multiplier_objective(u, &classes, pull);
                let (up, up_slope, _) = log_multiplier_objective(u + step, &classes, pull);
                let (down, down_slope, _) = log_multiplier_objective(u - step, &classes, pull);
                let numeric_slope = (up - down) / (2.0 * step);
                let numeric_curvature = (up_slope - down_slope) / (2.0 * step);
                assert!(
                    (slope - numeric_slope).abs() < 1e-5 * slope.abs().max(1.0),
                    "slope {slope} against {numeric_slope} at u {u}, pull {pull}"
                );
                assert!(
                    (curvature - numeric_curvature).abs() < 1e-5 * curvature.abs().max(1.0),
                    "curvature {curvature} against {numeric_curvature} at u {u}, pull {pull}"
                );
                assert!(curvature < 0.0, "concave in ln m");
            }
        }
    }

    /// **The pull's two limits**: a read group with no reads comes out at exactly one; one with a
    /// great many reads at the ratio its reads say, the pull's share of the evidence negligible.
    #[test]
    fn with_no_reads_the_pull_holds_one_and_with_many_the_reads_decide() {
        assert_eq!(fit_multiplier(&[], PULL_TOWARD_ONE_SLIPPED_READS, 1.7), 1.0);
        // Read at 2.5 times the class rate of 0.02 over a million reads: 50,000 slipped.
        let classes = [evidence(0.02, 950_000.0, 50_000.0)];
        let fitted = fit_multiplier(&classes, PULL_TOWARD_ONE_SLIPPED_READS, 1.0);
        assert!((fitted - 2.5).abs() < 2.5e-3, "fitted {fitted} against 2.5");
        // Without the pull, the reads' own ratio exactly.
        let alone = fit_multiplier(&classes, 0.0, 1.0);
        assert!((alone - 2.5).abs() < 1e-9, "fitted {alone} against 2.5");
    }

    /// **A multiplier never lifts a class past the ceiling**: a read group whose reads at the class
    /// of highest rate all slipped stops where that class reaches `MULTIPLIED_LEVEL_CEILING`.
    #[test]
    fn a_multiplier_stops_where_it_lifts_a_class_to_the_ceiling() {
        let classes = [evidence(0.5, 0.0, 100.0), evidence(0.01, 99.0, 1.0)];
        let fitted = fit_multiplier(&classes, 0.0, 1.0);
        assert!(fitted * 0.5 <= MULTIPLIED_LEVEL_CEILING, "fitted {fitted}");
        assert!(
            fitted * 0.5 > MULTIPLIED_LEVEL_CEILING * 0.999,
            "fitted {fitted}"
        );
    }

    /// **Scaled, the multipliers' weighted geometric mean is one** and every read group's slip
    /// rate is unchanged once the class rates take the factor.
    #[test]
    fn scaled_multipliers_have_a_weighted_geometric_mean_of_one() {
        let mut multipliers = vec![0.5, 2.0, 4.0, 1.0];
        let weights = [100.0, 300.0, 50.0, 0.0];
        let before = multipliers.clone();
        let factor = scale_to_a_typical_read_group(&mut multipliers, &weights);
        let mean_log: f64 = multipliers
            .iter()
            .zip(&weights)
            .map(|(m, w)| w * m.ln())
            .sum::<f64>()
            / weights.iter().sum::<f64>();
        assert!(mean_log.abs() < 1e-12, "weighted mean log {mean_log}");
        for (scaled, old) in multipliers.iter().zip(&before) {
            assert!((scaled * factor - old).abs() < 1e-12);
        }
    }

    /// **Against drawn truth, the multipliers land within two of their errors of it about 95 times
    /// in 100**, with the classes held at their true rates (spec §8.1). 200 read groups, their true
    /// multipliers log-normal with a standard deviation of 0.4 as in the stuttering simulation, at
    /// 3 reads a tract and at 30, 40 tracts a class over 20 classes. With no pull, so the errors are
    /// the reads' own; with the pull, the estimates move toward one and the 3-read case is checked
    /// only for that.
    #[test]
    fn drawn_multipliers_are_recovered_within_their_errors() {
        let mut draw = Draw(2026_11_02);
        for depth in [3, 30] {
            let read_groups = 200;
            let mut covered = 0;
            for _ in 0..read_groups {
                let truth = (0.4 * draw.normal()).exp();
                let (classes, by_tract) = drawn_read_group(truth, 40, depth, &mut draw);
                let fitted = fit_multiplier(&classes, 0.0, 1.0);
                let error = log_multiplier_error(fitted, &classes, &by_tract, 0.0)
                    .expect("every drawn read group has slipped reads");
                if (fitted.ln() - truth.ln()).abs() <= 2.0 * error {
                    covered += 1;
                }
            }
            assert!(
                (180..=200).contains(&covered),
                "{covered} of {read_groups} read groups within two errors at {depth} reads a tract"
            );
        }
        // With the pull, a read group of few slipped reads moves toward one, never past it.
        let truth = 3.0;
        let (classes, _) = drawn_read_group(truth, 2, 3, &mut draw);
        let alone = fit_multiplier(&classes, 0.0, 1.0);
        let pulled = fit_multiplier(&classes, PULL_TOWARD_ONE_SLIPPED_READS, 1.0);
        assert!(
            pulled < alone && pulled > 1.0,
            "pulled {pulled} between one and the reads' own {alone}"
        );
    }

    /// **Where tracts slip more or less than their class, the error still covers the truth**,
    /// because it is spread over tracts and not read off the curvature alone: each tract's rate its
    /// class's times a log-normal factor of log standard deviation 0.8, at 30 reads a tract. An
    /// error from the curvature alone, which assumes every read independent, held the truth for 154
    /// of these 200 read groups (measured 2026-11-02).
    #[test]
    fn the_error_allows_for_tracts_that_slip_more_than_their_class() {
        let mut draw = Draw(7);
        let read_groups = 200;
        let mut covered = 0;
        for _ in 0..read_groups {
            let truth = (0.4 * draw.normal()).exp();
            let (classes, by_tract) =
                drawn_read_group_with_tracts_apart(truth, 40, 30, 0.8, &mut draw);
            let fitted = fit_multiplier(&classes, 0.0, 1.0);
            let error = log_multiplier_error(fitted, &classes, &by_tract, 0.0)
                .expect("every drawn read group has slipped reads");
            if (fitted.ln() - truth.ln()).abs() <= 2.0 * error {
                covered += 1;
            }
        }
        assert!(
            (180..=200).contains(&covered),
            "{covered} of {read_groups} read groups within two errors"
        );
    }

    /// **No error for a read group with fewer slipped reads than the pull's pseudo-reads.**
    #[test]
    fn a_multiplier_resting_on_the_pull_carries_no_error() {
        let classes = [evidence(0.01, 1_000.0, 10.0)];
        let by_tract = vec![classes.to_vec()];
        assert_eq!(
            log_multiplier_error(1.0, &classes, &by_tract, PULL_TOWARD_ONE_SLIPPED_READS),
            None
        );
        assert!(log_multiplier_error(1.0, &classes, &by_tract, 5.0).is_some());
    }
}
