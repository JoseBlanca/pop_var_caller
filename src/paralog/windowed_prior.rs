//! **A duplication rate for each stretch of the genome, instead of one for the whole run.**
//! Research only (branch `research/giab-alpha`); switched on by `--hidden-dup-window`, with two
//! knobs in environment switches ([`crate::exp::windowed_paralog_prior`]). At a window of `0`,
//! nothing here runs.
//!
//! **What it is for.** The filter turns each record's likelihood ratio into a probability by
//! adding the log-odds of one run-wide rate, π, the share of records that are hidden
//! duplications ([`super::prior`]). Duplications cluster, so one rate is wrong almost everywhere:
//! on GIAB HG002 chr22 the same estimate gives about 0.15 over all records, about 0.0005 inside
//! GIAB's trusted regions and about 0.4 outside them. In easy ground the run-wide rate inflates
//! every probability, and the filter removes real variants there.
//!
//! **What this does instead.** The contig is cut into bins of 10,000 bases. For each bin that holds
//! a scored record, the rate is re-estimated from the records within half a window of the bin's
//! centre, by the same iteration as the run-wide rate, but pulled towards the run-wide rate as if
//! `K` extra records sat in the window at that rate. A window with few records therefore keeps
//! about the run-wide rate; a window with many keeps its own. In symbols, the bin's rate `π_b` is
//! the fixed point of
//!
//! ```text
//! π_b = ( Σ_i σ(LR_i + logit π_b) + K·π ) / (n + K)
//! ```
//!
//! where the sum runs over the `n` scored records within `W/2` of the bin's centre, `σ` is the
//! logistic function, `π` the run-wide rate and `K` the pull towards it, counted in records. Each
//! record's probability of being a duplication is then `σ(LR + logit π_b)` with its own bin's rate.
//!
//! **The cut** is one of two rules ([`WindowedCut`]): the run's tail false-discovery rule applied
//! to these local probabilities, or a fixed probability taken from today's cut.
//!
//! **Records with no finite ratio** — unscored ones — stay as they are today: they are not counted
//! in any window, carry no probability, and are never flagged.
//!
//! **Cost.** Pass two keeps a contig and a position (12 bytes) beside each record's ratio, and the
//! calibration briefly holds an index (8 bytes) a record to sort them, then one probability (8
//! bytes) a record until pass three. All of it is flat in the sample count. The time is one
//! iteration per occupied bin over the records in its window, so it grows with the window: about
//! `records × W / 10,000 × iterations` logistic evaluations, at most 50 iterations a bin.

use rayon::prelude::*;

use crate::float;

/// The width of a bin: every record in one bin shares one local rate.
pub const LOCAL_PRIOR_BIN_BASES: u64 = 10_000;
/// The iteration stops when the local rate moves by less than this.
pub const LOCAL_PRIOR_TOLERANCE: f64 = 1e-7;
/// The iteration's cap.
pub const LOCAL_PRIOR_MAX_ITERATIONS: usize = 50;
/// The local rate is held inside `[floor, 1 − floor]`, so its log-odds stay finite.
pub const LOCAL_PRIOR_FLOOR: f64 = 1e-9;

/// **How the windowed probabilities become a verdict.**
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WindowedCut {
    /// Rank the records by local probability, highest first, and flag the longest run from the
    /// top whose mean `1 − probability` is at or below the target false-discovery rate. With
    /// one rate everywhere this is today's tail-FDR cut, without the histogram's binning.
    TailFdr,
    /// Flag a record whose local probability exceeds the probability today's cut stands for:
    /// `σ(lr_cut + logit π)`, with today's cut and the run-wide rate.
    Posterior,
}

impl WindowedCut {
    /// The name the header and the switch use: `fdr` or `posterior`.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::TailFdr => "fdr",
            Self::Posterior => "posterior",
        }
    }
}

/// **The three research knobs, together.**
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowedPriorConfig {
    /// The window's width in bases: the records within half of it of a bin's centre set that
    /// bin's rate.
    pub window_bases: u64,
    /// How hard the local rate is pulled towards the run-wide one, in records.
    pub shrink_records: f64,
    /// Which rule turns the local probabilities into verdicts.
    pub cut: WindowedCut,
}

/// The lowest, median and highest local rate over the occupied bins.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LocalRateSpread {
    /// The lowest.
    pub lowest: f64,
    /// The median (the lower middle value on an even count).
    pub median: f64,
    /// The highest.
    pub highest: f64,
}

/// **What the windowed calibration settles for pass three.**
#[derive(Debug, Clone)]
pub struct WindowedVerdicts {
    /// The knobs it ran with.
    pub config: WindowedPriorConfig,
    /// One probability a record, **in spill order**, `NaN` where the record was not scored. It
    /// replaces `σ(LR + logit π)` as the record's `PARALOG_POST`.
    pub posteriors: Vec<f64>,
    /// The probability the cut resolved to, or `None` where nothing is flagged.
    pub posterior_cut: Option<f64>,
    /// How many bins held a scored record and so got a rate of their own.
    pub occupied_bins: usize,
    /// The spread of those rates, or `None` where no bin was occupied.
    pub local_rates: Option<LocalRateSpread>,
}

impl WindowedVerdicts {
    /// Whether a record with this local probability is flagged. A `NaN` — an unscored record —
    /// never is.
    ///
    /// The two rules differ at the boundary on purpose. The tail-FDR rule takes whole groups of
    /// tied probabilities and its cut is the lowest probability it took, so the cut itself is
    /// flagged; the probability rule flags what *exceeds* the cut.
    #[must_use]
    pub fn flags(&self, posterior: f64) -> bool {
        match (self.posterior_cut, self.config.cut) {
            (None, _) => false,
            (Some(cut), WindowedCut::TailFdr) => posterior >= cut,
            (Some(cut), WindowedCut::Posterior) => posterior > cut,
        }
    }
}

/// **One bin's rate**: the fixed point of `π_b = (Σ σ(LR_i + logit π_b) + K·π) / (n + K)`.
///
/// `ratios` are the finite ratios in the bin's window, `global_rate` is `π` and `shrink_records` is
/// `K`. The iteration starts at `π`, stops when the rate moves by less than
/// [`LOCAL_PRIOR_TOLERANCE`] or after [`LOCAL_PRIOR_MAX_ITERATIONS`], and holds the rate inside
/// `[LOCAL_PRIOR_FLOOR, 1 − LOCAL_PRIOR_FLOOR]`. A window with no records and no pull (`n + K = 0`)
/// has nothing to estimate from and keeps `π`.
#[must_use]
pub fn local_rate(ratios: &[f64], global_rate: f64, shrink_records: f64) -> f64 {
    let global_rate = clamp_rate(global_rate);
    let denominator = ratios.len() as f64 + shrink_records;
    if denominator <= 0.0 {
        return global_rate;
    }
    let pull = shrink_records * global_rate;
    let mut rate = global_rate;
    for _ in 0..LOCAL_PRIOR_MAX_ITERATIONS {
        let offset = logit(rate);
        let expected: f64 = ratios.iter().map(|&ratio| sigmoid(ratio + offset)).sum();
        let next = clamp_rate((expected + pull) / denominator);
        let settled = (next - rate).abs() < LOCAL_PRIOR_TOLERANCE;
        rate = next;
        if settled {
            break;
        }
    }
    rate
}

/// The probability `σ(LR + logit rate)` that a record with this ratio is a duplication.
#[must_use]
pub fn posterior_at(ratio: f64, rate: f64) -> f64 {
    sigmoid(ratio + logit(clamp_rate(rate)))
}

/// **The tail-FDR rule on a set of probabilities**: the lowest probability of the longest run,
/// from the highest probability down, whose mean `1 − probability` is at or below `target_fdr`.
/// `None` where not even the highest group meets it.
///
/// Tied probabilities are taken or left as one group, because a cut cannot separate them. Sorted
/// highest first, each added `1 − probability` is at least the mean before it, so the mean never
/// falls; the first group that lifts it past the target ends the run. Non-finite values are
/// ignored.
#[must_use]
pub fn tail_fdr_posterior_cut(posteriors: &[f64], target_fdr: f64) -> Option<f64> {
    let mut sorted: Vec<f64> = posteriors
        .iter()
        .copied()
        .filter(|p| p.is_finite())
        .collect();
    sorted.sort_unstable_by(|a, b| b.total_cmp(a));
    let mut cut = None;
    let mut taken = 0usize;
    let mut false_share_sum = 0.0;
    let mut start = 0;
    while start < sorted.len() {
        let value = sorted[start];
        let end = start + sorted[start..].iter().take_while(|&&p| p == value).count();
        let group_sum = (end - start) as f64 * (1.0 - value);
        let mean = (false_share_sum + group_sum) / (taken + end - start) as f64;
        if mean > target_fdr {
            break;
        }
        false_share_sum += group_sum;
        taken = end;
        cut = Some(value);
        start = end;
    }
    cut
}

/// **Every record's local probability, and the cut.**
///
/// `ratios`, `contigs` and `positions` are parallel, in spill order. `global_rate` is the run's
/// fitted rate (or its fallback), `global_lr_cut` the ratio today's cut resolved to, and
/// `target_fdr` the operator's target.
///
/// # Panics
///
/// If the three slices differ in length — pass two fills them together, one entry a record.
#[must_use]
pub fn windowed_verdicts(
    ratios: &[f64],
    contigs: &[u32],
    positions: &[u64],
    global_rate: f64,
    global_lr_cut: Option<f64>,
    target_fdr: f64,
    config: &WindowedPriorConfig,
) -> WindowedVerdicts {
    assert!(
        ratios.len() == contigs.len() && ratios.len() == positions.len(),
        "one contig and one position a ratio"
    );

    // The scored records, by contig and then position. The spill is already in that order; the
    // sort is what lets this not depend on it.
    let mut order: Vec<usize> = (0..ratios.len())
        .filter(|&i| ratios[i].is_finite())
        .collect();
    order.sort_by_key(|&i| (contigs[i], positions[i]));

    let mut posteriors = vec![f64::NAN; ratios.len()];
    let mut all_rates: Vec<f64> = Vec::new();
    let half_window = config.window_bases / 2;

    let mut contig_start = 0;
    while contig_start < order.len() {
        let contig = contigs[order[contig_start]];
        let contig_end = contig_start
            + order[contig_start..]
                .iter()
                .take_while(|&&i| contigs[i] == contig)
                .count();
        let members = &order[contig_start..contig_end];
        let contig_positions: Vec<u64> = members.iter().map(|&i| positions[i]).collect();
        let contig_ratios: Vec<f64> = members.iter().map(|&i| ratios[i]).collect();

        let mut bins: Vec<u64> = contig_positions
            .iter()
            .map(|&p| p / LOCAL_PRIOR_BIN_BASES)
            .collect();
        bins.dedup();

        let rates: Vec<f64> = bins
            .par_iter()
            .map(|&bin| {
                let centre = bin * LOCAL_PRIOR_BIN_BASES + LOCAL_PRIOR_BIN_BASES / 2;
                let low = centre.saturating_sub(half_window);
                let high = centre.saturating_add(half_window);
                let first = contig_positions.partition_point(|&p| p < low);
                let past = contig_positions.partition_point(|&p| p <= high);
                local_rate(
                    &contig_ratios[first..past],
                    global_rate,
                    config.shrink_records,
                )
            })
            .collect();

        // Bins and records are both in position order, so one walk pairs each record with its
        // bin's rate.
        let mut bin_at = 0;
        for (k, &record) in members.iter().enumerate() {
            let bin = contig_positions[k] / LOCAL_PRIOR_BIN_BASES;
            while bins[bin_at] != bin {
                bin_at += 1;
            }
            posteriors[record] = posterior_at(contig_ratios[k], rates[bin_at]);
        }
        all_rates.extend(rates);
        contig_start = contig_end;
    }

    let posterior_cut = match config.cut {
        WindowedCut::TailFdr => tail_fdr_posterior_cut(&posteriors, target_fdr),
        WindowedCut::Posterior => global_lr_cut.map(|lr| posterior_at(lr, global_rate)),
    };

    let occupied_bins = all_rates.len();
    all_rates.sort_unstable_by(f64::total_cmp);
    let local_rates = (!all_rates.is_empty()).then(|| LocalRateSpread {
        lowest: all_rates[0],
        median: all_rates[(all_rates.len() - 1) / 2],
        highest: all_rates[all_rates.len() - 1],
    });

    WindowedVerdicts {
        config: *config,
        posteriors,
        posterior_cut,
        occupied_bins,
        local_rates,
    }
}

fn clamp_rate(rate: f64) -> f64 {
    rate.clamp(LOCAL_PRIOR_FLOOR, 1.0 - LOCAL_PRIOR_FLOOR)
}

fn logit(p: f64) -> f64 {
    float::ln(p / (1.0 - p))
}

fn sigmoid(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + float::exp(-x))
    } else {
        let e = float::exp(x);
        e / (1.0 + e)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::paralog::{
        CalibrationConfig, ParalogLrHistogram, calibrate_from_the_ratio_histogram,
    };

    /// The fixed point solved independently, by bisection on `f(p) = rhs(p) − p`, which is
    /// positive below the fixed point and negative above it.
    fn fixed_point_by_bisection(ratios: &[f64], global: f64, shrink: f64) -> f64 {
        let rhs = |p: f64| {
            let s: f64 = ratios.iter().map(|&r| sigmoid(r + logit(p))).sum();
            (s + shrink * global) / (ratios.len() as f64 + shrink)
        };
        let (mut low, mut high) = (1e-9, 1.0 - 1e-9);
        for _ in 0..200 {
            let mid = 0.5 * (low + high);
            if rhs(mid) > mid {
                low = mid;
            } else {
                high = mid;
            }
        }
        0.5 * (low + high)
    }

    /// **The local rate on four hand-made records.** Two sit at a ratio of +10 and two at −10.
    /// With no pull the fixed point is exactly one half, since `σ(x) + σ(−x) = 1`. With a pull of
    /// four records at a run-wide 0.1 it lies between the two, where an independent bisection puts
    /// it. With no records the pull is all there is, and the rate is the run-wide one.
    #[test]
    fn the_local_rate_is_the_fixed_point_of_the_shrunk_iteration() {
        let ratios = [10.0, 10.0, -10.0, -10.0];
        let unpulled = local_rate(&ratios, 0.1, 0.0);
        assert!((unpulled - 0.5).abs() < 1e-6, "no pull: {unpulled}");

        let pulled = local_rate(&ratios, 0.1, 4.0);
        let expected = fixed_point_by_bisection(&ratios, 0.1, 4.0);
        assert!(pulled > 0.1 && pulled < 0.5, "pulled: {pulled}");
        assert!(
            (pulled - expected).abs() < 1e-5,
            "pulled {pulled}, bisection {expected}"
        );

        assert!((local_rate(&[], 0.1, 10.0) - 0.1).abs() < 1e-12);
        assert!((local_rate(&[], 0.1, 0.0) - 0.1).abs() < 1e-12);
    }

    /// **With one rate everywhere, the tail-FDR rule on probabilities flags what today's
    /// histogram cut flags**, apart from records whose ratio lies within a bin or so of the cut
    /// — the histogram decides by bin, this by record.
    #[test]
    fn the_tail_fdr_cut_on_one_rate_is_todays_cut() {
        // 9,000 records spread over a real-variant range and 1,000 over a duplication range, with
        // a stretch between so the cut falls among records.
        let mut ratios: Vec<f64> = Vec::new();
        for i in 0..9_000 {
            ratios.push(-15.0 + 0.002 * i as f64); // [-15, 3)
        }
        for i in 0..1_000 {
            ratios.push(0.0 + 0.012 * i as f64); // [0, 12)
        }
        let mut histogram = ParalogLrHistogram::with_defaults();
        for &r in &ratios {
            histogram.push(r);
        }
        for target in [0.01, 0.05, 0.2] {
            let today = calibrate_from_the_ratio_histogram(
                &histogram,
                target,
                &CalibrationConfig::default(),
            );
            let rate = today.prior.prior_probability;
            let lr_cut = today.lr_threshold.expect("the target is reachable");
            let posteriors: Vec<f64> = ratios.iter().map(|&r| posterior_at(r, rate)).collect();
            let cut = tail_fdr_posterior_cut(&posteriors, target).expect("reachable here too");

            let mut disagree = 0;
            for (&r, &p) in ratios.iter().zip(&posteriors) {
                if today.flags(r) != (p >= cut) {
                    disagree += 1;
                    assert!(
                        (r - lr_cut).abs() < 0.2,
                        "target {target}: ratio {r} is far from today's cut {lr_cut} and the two \
                         rules disagree on it"
                    );
                }
            }
            assert!(
                disagree < ratios.len() / 100,
                "target {target}: {disagree} records disagree"
            );
        }
    }

    /// Tied probabilities are taken whole or not at all, and a target nothing meets flags nothing.
    #[test]
    fn the_tail_fdr_cut_takes_ties_whole() {
        // Highest first: 0.99 alone has a mean false share of 0.01; with both 0.97s it is
        // (0.01 + 0.03 + 0.03) / 3 ≈ 0.0233; with 0.5 too, 0.1425. Taking one 0.97 of the two
        // would give 0.02, which a target of 0.021 admits, and the rule must not do that.
        let posteriors = [0.99, 0.97, 0.97, 0.5, f64::NAN];
        assert_eq!(tail_fdr_posterior_cut(&posteriors, 0.011), Some(0.99));
        assert_eq!(tail_fdr_posterior_cut(&posteriors, 0.021), Some(0.99));
        assert_eq!(tail_fdr_posterior_cut(&posteriors, 0.025), Some(0.97));
        assert_eq!(tail_fdr_posterior_cut(&posteriors, 0.2), Some(0.5));
        assert_eq!(tail_fdr_posterior_cut(&posteriors, 0.005), None);
    }

    /// **A record among many duplication-shaped neighbours is more likely a duplication than an
    /// identical record alone.** Twenty records with a ratio of 5 sit within 10 kb of each other;
    /// one more with the same ratio sits 5 Mb away; ten thousand real-variant-shaped records fill
    /// the contig around both.
    #[test]
    fn a_dense_window_raises_its_records_probability() {
        let mut ratios = Vec::new();
        let mut positions = Vec::new();
        for i in 0..10_000u64 {
            ratios.push(-10.0);
            positions.push(i * 1_000 + 500);
        }
        for i in 0..20u64 {
            ratios.push(5.0);
            positions.push(100_000 + i * 400);
        }
        ratios.push(5.0);
        positions.push(5_000_123);
        let contigs = vec![0u32; ratios.len()];

        let mut histogram = ParalogLrHistogram::with_defaults();
        for &r in &ratios {
            histogram.push(r);
        }
        let today =
            calibrate_from_the_ratio_histogram(&histogram, 0.01, &CalibrationConfig::default());
        let rate = today.prior.prior_probability;
        let config = WindowedPriorConfig {
            window_bases: 50_000,
            shrink_records: 10.0,
            cut: WindowedCut::TailFdr,
        };
        let windowed = windowed_verdicts(
            &ratios,
            &contigs,
            &positions,
            rate,
            today.lr_threshold,
            0.01,
            &config,
        );
        let dense = windowed.posteriors[10_000];
        let isolated = windowed.posteriors[10_020];
        assert!(
            dense > isolated,
            "dense {dense} should exceed isolated {isolated}"
        );
        assert!(dense > posterior_at(5.0, rate), "dense {dense} vs run-wide");
        assert_eq!(windowed.occupied_bins, 1_000);
    }

    /// An unscored record carries no probability and is never flagged, and it does not count in
    /// its window.
    #[test]
    fn an_unscored_record_has_no_probability_and_is_not_counted() {
        let config = WindowedPriorConfig {
            window_bases: 50_000,
            shrink_records: 0.0,
            cut: WindowedCut::TailFdr,
        };
        let with = windowed_verdicts(
            &[f64::NAN, 10.0, -10.0],
            &[0, 0, 0],
            &[10, 20, 30],
            0.1,
            None,
            0.5,
            &config,
        );
        let without =
            windowed_verdicts(&[10.0, -10.0], &[0, 0], &[20, 30], 0.1, None, 0.5, &config);
        assert!(with.posteriors[0].is_nan());
        assert!(!with.flags(with.posteriors[0]));
        assert_eq!(with.posteriors[1], without.posteriors[0]);
        assert_eq!(with.posteriors[2], without.posteriors[1]);
    }

    /// The probability rule flags what exceeds the probability today's cut stands for, and with
    /// no cut flags nothing.
    #[test]
    fn the_posterior_rule_uses_todays_cut_at_the_run_wide_rate() {
        let config = WindowedPriorConfig {
            window_bases: 50_000,
            shrink_records: 10.0,
            cut: WindowedCut::Posterior,
        };
        let windowed = windowed_verdicts(
            &[3.0, -3.0],
            &[0, 1],
            &[5, 5],
            0.2,
            Some(1.5),
            0.01,
            &config,
        );
        assert_eq!(windowed.posterior_cut, Some(posterior_at(1.5, 0.2)));
        let none = windowed_verdicts(&[3.0], &[0], &[5], 0.2, None, 0.01, &config);
        assert_eq!(none.posterior_cut, None);
        assert!(!none.flags(1.0));
    }
}
