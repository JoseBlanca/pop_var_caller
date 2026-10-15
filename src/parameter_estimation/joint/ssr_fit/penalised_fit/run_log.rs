//! **What the repeat-tract fit prints as it runs** (spec `str_penalised_fit.md` §2.9), so that a
//! person reading the log of `estimate-parameters` can tell how each motif length's fit went
//! without opening the parameters file.
//!
//! While a motif length is being fitted, **at most once a minute, where it has got to**: the round
//! reached and the largest change any number made in it, as every looping stage of
//! `estimate-parameters` does.
//!
//! After each motif length is fitted:
//!
//! - **one line on how the fit ran**: how many rounds, whether it settled or used them all, how F
//!   was read, and the time in each step: weighing genotypes and counting (steps 1 and 2 of each
//!   round), drawing the curves (step 3, with the part of it spent in the rounds that also chose
//!   how stiff the curves are), reading F off pairs of samples, and the standard errors. It is
//!   printed as a progress line, so it also carries the time into the stage and, on Linux, the
//!   process's memory now and the most it has held so far (`VmHWM`). On macOS nothing is printed
//!   for memory. F is named there, the first time, as "the between-locus share": how much the
//!   allele-length frequencies vary from one locus to the next;
//! - **how stiff each of the five kinds of curve came out**: slip rate, shorter share, fall-off, F
//!   and the allele-length shares (one strength for all of them), against the range searched. With
//!   one or two repeat counts there is no bend to penalise, and the line says so instead;
//! - **which strengths could not be scored**, printed only when some could not: a strength whose
//!   curve's fit could not be factorised is never chosen, and a search that scores none keeps the
//!   strength held before;
//! - **where bends were eased**: the repeat counts where a curve's reads insisted on a kink, so its
//!   penalty was eased there;
//! - **why standard errors are absent**, printed only when they are for this reason: the curvature
//!   they are read off could not be inverted, so they are written as not placed;
//! - **one line a class**: its loci with reads, its (locus, sample) pairs with reads and how many
//!   distinct rows of read counts they hold, and whether F is absent and why.
//!
//! At the end, **how many standard errors were withheld** because too few loci stand behind them
//! (spec §2.6, §7). With no class at all, one line says there was nothing to fit.
//!
//! Every line is built here from what the fit already returned; nothing here changes a number.

use std::collections::BTreeSet;
use std::time::Duration;

use super::*;

/// One line of the repeat-tract fit's log, and how it is to be printed.
#[derive(Debug, Clone, PartialEq)]
pub(super) enum Line {
    /// A line that closes a step, printed as a progress line: with the time into the stage and,
    /// on Linux, the process's memory.
    StepEnded(String),
    /// A detail of the step just closed, printed as it is.
    Detail(String),
}

impl Line {
    /// The line's text, however it is printed.
    pub(super) fn text(&self) -> &str {
        match self {
            Self::StepEnded(text) | Self::Detail(text) => text,
        }
    }
}

/// `n` followed by the noun for one or for many: `1 class`, `5 classes`.
fn count_of(n: usize, one: &str, many: &str) -> String {
    format!("{n} {}", if n == 1 { one } else { many })
}

/// The one line printed when no repeat class is handed in.
pub(super) fn no_classes_line() -> String {
    "no repeat classes were handed in, so there was nothing to fit".to_owned()
}

/// The line that opens the stage: how many classes, of how many motif lengths.
pub(super) fn opening_line(strata: &[StratumEvidence]) -> String {
    let motif_lengths: BTreeSet<u8> = strata
        .iter()
        .map(|evidence| evidence.stratum.period)
        .collect();
    format!(
        "fitting {classes} (one motif length and one reference repeat count each) of {motif_lengths} \
         by curves along the repeat count, one motif length at a time",
        classes = count_of(strata.len(), "repeat class", "repeat classes"),
        motif_lengths = count_of(motif_lengths.len(), "motif length", "motif lengths"),
    )
}

/// **Where a motif length's fit has got to**, for the line printed at most once a minute while it
/// runs: the round just ended and the largest change any number made in it.
pub(super) fn round_line(
    period: u8,
    round: usize,
    largest_change: f64,
    config: &PenalisedFitConfig,
) -> String {
    format!(
        "motif length {period}, round {round} of at most {max_rounds}: the largest change in any \
         number was {largest_change:.1e} (it settles once every change is below {tolerance:e})",
        max_rounds = config.max_rounds,
        tolerance = config.tolerance,
    )
}

/// **Every line one motif length's fit gives the log**, in the order they are printed.
///
/// `classes` are the motif length's classes in the order they were handed to its fit, and
/// `outcomes` and `distinct_rows` are theirs in the same order; `fit` is what the motif length's
/// fit returned, its classes already taken out.
pub(super) fn lines_of_a_motif_length(
    classes: &[&StratumEvidence],
    outcomes: &[StratumOutcome],
    distinct_rows: &[Option<usize>],
    fit: &PenalisedPeriodFit,
    config: &PenalisedFitConfig,
) -> Vec<Line> {
    let Some(first) = classes.first() else {
        return Vec::new();
    };
    let period = first.stratum.period;
    let Some(curves) = fit.curves.as_ref() else {
        return vec![Line::StepEnded(format!(
            "motif length {period}: no reads at any locus in its {classes}, so nothing was fitted",
            classes = count_of(classes.len(), "class", "classes"),
        ))];
    };
    let fitted = outcomes
        .iter()
        .filter(|outcome| matches!(outcome, StratumOutcome::Fitted(_)))
        .count();
    let mut lines = vec![Line::StepEnded(how_the_fit_ran(
        period, outcomes, fitted, fit, curves, config,
    ))];
    // With one or two repeat counts with reads the curve through them is a line whatever the
    // strength, so the strength held is the top of the range searched and means nothing (spec
    // §2.3). Every fitted class has reads, so it is a repeat count with reads.
    if fitted <= 2 {
        lines.push(Line::Detail(format!(
            "motif length {period}, how stiff each curve came out: {}, so no curve has a bend to \
             penalise; each class keeps its own counts",
            if fitted == 1 {
                "one repeat count with reads"
            } else {
                "only two repeat counts with reads"
            },
        )));
    } else {
        lines.push(Line::Detail(stiffness_line(period, curves)));
    }
    // A fit kept at a strength it could not factorise can happen with one or two repeat counts
    // too, so this line is not tied to a curve having bends.
    lines.extend(unscored_strengths_line(period, curves).map(Line::Detail));
    if fitted > 2 {
        lines.push(Line::Detail(eased_bends_line(
            period,
            curves,
            config.allele_span,
        )));
    }
    if let Some(undetermined) = fit.undetermined_errors {
        lines.push(Line::Detail(undetermined_errors_line(period, undetermined)));
    }
    let f_measurable = curves.between_locus_share.is_some();
    lines.extend(classes.iter().zip(outcomes).zip(distinct_rows).map(
        |((evidence, outcome), distinct)| {
            Line::Detail(class_line(evidence, outcome, *distinct, f_measurable))
        },
    ));
    lines
}

/// The motif length's rounds, how it stopped, how F was read and the time in each step.
fn how_the_fit_ran(
    period: u8,
    outcomes: &[StratumOutcome],
    fitted: usize,
    fit: &PenalisedPeriodFit,
    curves: &SmoothCurves,
    config: &PenalisedFitConfig,
) -> String {
    // A reading of F after the rounds ran out is followed by no round, so the last change counted
    // may predate F's last move: only a fit that stopped short of its limit is called settled.
    // F's change is counted at a hundredth of its logarithm's move, hence "about".
    let ending = if fit.rounds < config.max_rounds && fit.last_change < config.tolerance {
        format!(
            "settled (no number moved by more than about {:e} in the last round)",
            config.tolerance
        )
    } else {
        format!(
            "used all {} rounds; the largest change in the last one was {:.1e}",
            config.max_rounds, fit.last_change
        )
    };
    let f = "F (the between-locus share: how much the allele-length frequencies vary from one \
             locus to the next)";
    let how_f_was_read = if curves.between_locus_share.is_none() {
        format!("{f} absent: no locus has two samples with reads")
    } else if fit.pair_alternations > 0 {
        format!(
            "{f} read off pairs of samples at a locus {}",
            count_of(fit.pair_alternations, "time", "times")
        )
    } else {
        format!("{f} read off each sample's own two copies only")
    };
    let time = |seconds: f64| progress::duration(Duration::from_secs_f64(seconds.max(0.0)));
    let which = if outcomes.len() == 1 {
        "its one class fitted".to_owned()
    } else {
        format!("{fitted} of {} classes fitted together", outcomes.len())
    };
    format!(
        "motif length {period}: {which} in {rounds}, {ending}; {how_f_was_read}; time weighing \
         genotypes and counting {expectation}, drawing the curves {maximisation} \
         ({strength_search} of it in the {choosing} that also chose how stiff they are), reading \
         F off pairs {pairs}, standard errors {errors}",
        rounds = count_of(fit.rounds, "round", "rounds"),
        expectation = time(fit.expectation_seconds),
        maximisation = time(fit.maximisation_seconds),
        strength_search = time(fit.strength_search_seconds),
        choosing = count_of(fit.rounds_choosing_strength, "round", "rounds"),
        pairs = time(fit.pair_seconds),
        errors = time(fit.error_seconds),
    )
}

/// A penalty strength as a reader takes it in: two significant figures, written plainly from
/// 0.01 to below 1,000 and as a power of ten outside that.
///
/// Built from the digits of the rounded scientific form (`2.4e0`), so the plain form shows
/// exactly the two figures rounded, with no float arithmetic of its own.
fn readable(strength: f64) -> String {
    let scientific = format!("{strength:.1e}");
    let parsed = scientific.split_once('e').and_then(|(mantissa, exponent)| {
        let (whole, tenth) = mantissa.split_once('.')?;
        Some((
            whole.to_owned(),
            tenth.to_owned(),
            exponent.parse::<i32>().ok()?,
        ))
    });
    let Some((whole, tenth, exponent)) = parsed else {
        return scientific;
    };
    match exponent {
        0 => format!("{whole}.{tenth}"),
        1 | 2 => format!("{whole}{tenth}{}", "0".repeat(exponent as usize - 1)),
        -2 | -1 => format!(
            "0.{}{whole}{tenth}",
            "0".repeat(exponent.unsigned_abs() as usize - 1)
        ),
        _ => scientific,
    }
}

/// How stiff each of the five kinds of curve came out: its penalty strength on bends, against the
/// range searched.
fn stiffness_line(period: u8, curves: &SmoothCurves) -> String {
    let strength = |curve: &ParameterCurve| readable(curve.smoothing);
    format!(
        "motif length {period}, how stiff each curve came out, as the penalty on its bends \
         (searched from {loosest:.0e}, where each repeat count stands on its own, to \
         {stiffest:.0e}, a straight line): slip rate {level}, shorter share {shorter}, fall-off \
         {fall_off}, F {f}, allele-length shares (one strength for all of them) {alleles}",
        loosest = float::exp(LN_SMOOTHING_RANGE.0),
        stiffest = float::exp(LN_SMOOTHING_RANGE.1),
        level = strength(&curves.level),
        shorter = strength(&curves.shorter_share),
        fall_off = strength(&curves.fall_off),
        f = curves
            .between_locus_share
            .as_ref()
            .map_or_else(|| "absent".to_owned(), strength),
        alleles = curves
            .allele_length_variation
            .first()
            .map_or_else(|| "none".to_owned(), strength),
    )
}

/// **Which curves met strengths that could not be scored**, for the line printed only when one
/// did: a strength whose fit could not be factorised, even with a larger ridge, is never chosen;
/// a search that could score none keeps the strength held before; and a fit at a strength no
/// search chose that could not be factorised keeps the values Newton's method had reached. `None`
/// when every strength tried was scored and every fit kept was factorised.
fn unscored_strengths_line(period: u8, curves: &SmoothCurves) -> Option<String> {
    let mut named: Vec<(&str, &ParameterCurve)> = vec![
        ("slip rate", &curves.level),
        ("shorter share", &curves.shorter_share),
        ("fall-off", &curves.fall_off),
    ];
    if let Some(curve) = &curves.between_locus_share {
        named.push(("F", curve));
    }
    if let Some(curve) = curves.allele_length_variation.first() {
        named.push(("allele-length shares", curve));
    }
    let rejected: Vec<String> = named
        .iter()
        .filter(|(_, curve)| curve.strengths_rejected > 0)
        .map(|(name, curve)| format!("{name} {}", curve.strengths_rejected))
        .collect();
    let kept: Vec<String> = named
        .iter()
        .filter(|(_, curve)| curve.searches_without_a_strength > 0)
        .map(|(name, curve)| {
            format!(
                "{name} in {}",
                count_of(curve.searches_without_a_strength, "search", "searches")
            )
        })
        .collect();
    let unfactorised: Vec<String> = named
        .iter()
        .filter(|(_, curve)| curve.unfactorised_fits_kept > 0)
        .map(|(name, curve)| format!("{name} {}", curve.unfactorised_fits_kept))
        .collect();
    if rejected.is_empty() && kept.is_empty() && unfactorised.is_empty() {
        return None;
    }
    let mut parts = Vec::new();
    if !rejected.is_empty() {
        parts.push(format!(
            "strengths never chosen because the curve's fit there could not be factorised, even \
             with a larger ridge: {}",
            rejected.join(", ")
        ));
    }
    if !kept.is_empty() {
        parts.push(format!(
            "no strength tried could be scored, so the strength held before was kept: {}",
            kept.join(", ")
        ));
    }
    if !unfactorised.is_empty() {
        parts.push(format!(
            "fits at a strength no search chose that could not be factorised, so the curve kept \
             the values Newton's method had reached when it stopped: {}",
            unfactorised.join(", ")
        ));
    }
    Some(format!(
        "motif length {period}, strengths that could not be scored: {}",
        parts.join("; ")
    ))
}

/// **Why some of a motif length's standard errors are absent**: the curvature they are read off
/// could not be inverted ([`UndeterminedErrors`]). Every error so left out, and not already
/// withheld for too few loci, is written as not placed.
fn undetermined_errors_line(period: u8, undetermined: UndeterminedErrors) -> String {
    let (which, curvature) = match undetermined {
        UndeterminedErrors::EveryNumber => (
            "none determined for any of its numbers",
            "the fit's curvature, with the curves' bend penalties added",
        ),
        UndeterminedErrors::EveryNumberButBetweenLocusShare => (
            "none determined for any of its numbers except F, which is read off pairs of samples \
             with a curvature of its own",
            "the fit's curvature, with the curves' bend penalties added",
        ),
        UndeterminedErrors::BetweenLocusShare => (
            "none determined for F",
            "the curvature of F read off pairs of samples, with its curve's bend penalty added",
        ),
    };
    format!(
        "motif length {period}, standard errors: {which}: {curvature}, could not be inverted, \
         either because some combination of numbers has next to no reads behind it (the slope of \
         a stiff curve through a single repeat count's slipped reads, for one) or because the \
         curvature is not positive in some direction at the values the fit reached; those errors \
         are written as not placed, except where too few loci already withheld them"
    )
}

/// The curve of the share of alleles `offset` repeats from the reference length, in words.
fn allele_share_name(offset: i32) -> String {
    match offset {
        0 => "share of alleles at the reference length".to_owned(),
        _ => format!(
            "share of alleles {} {}",
            count_of(offset.unsigned_abs() as usize, "repeat", "repeats"),
            if offset > 0 { "longer" } else { "shorter" },
        ),
    }
}

/// Where each curve was let bend because its reads insisted on a kink. An allele-length share is
/// named by how many repeats its alleles lie from the reference length.
fn eased_bends_line(period: u8, curves: &SmoothCurves, allele_span: i32) -> String {
    let mut eased: Vec<(String, &[u64])> = vec![
        ("slip rate".to_owned(), curves.level.eased_at.as_slice()),
        (
            "shorter share".to_owned(),
            curves.shorter_share.eased_at.as_slice(),
        ),
        ("fall-off".to_owned(), curves.fall_off.eased_at.as_slice()),
    ];
    if let Some(curve) = &curves.between_locus_share {
        eased.push(("F".to_owned(), curve.eased_at.as_slice()));
    }
    for (allele, curve) in curves.allele_length_variation.iter().enumerate() {
        let offset = allele as i32 - allele_span;
        eased.push((allele_share_name(offset), curve.eased_at.as_slice()));
    }
    let named: Vec<String> = eased
        .into_iter()
        .filter(|(_, at)| !at.is_empty())
        .map(|(curve, at)| {
            let repeats: Vec<String> = at.iter().map(u64::to_string).collect();
            format!("{curve} at {} repeats", repeats.join(", "))
        })
        .collect();
    if named.is_empty() {
        format!("motif length {period}, bends eased: none; no curve's reads insisted on a kink")
    } else {
        format!(
            "motif length {period}, bends eased because the reads insisted on a kink: {}",
            named.join("; ")
        )
    }
}

/// One class's line: its reads, how many distinct rows they hold, and whether F is absent.
///
/// `f_measurable` says whether any locus of the motif length holds two samples with reads; when
/// none does, every class's F is absent for that reason.
fn class_line(
    evidence: &StratumEvidence,
    outcome: &StratumOutcome,
    distinct_rows: Option<usize>,
    f_measurable: bool,
) -> String {
    let Stratum {
        period,
        reference_repeats,
    } = evidence.stratum;
    let StratumOutcome::Fitted(fit) = outcome else {
        return format!(
            "motif length {period} at {reference_repeats} repeats: no reads at any locus, so not \
             fitted"
        );
    };
    let distinct = distinct_rows.map_or_else(
        || "distinct rows not counted (no round ran)".to_owned(),
        |distinct| {
            format!(
                "holding {} of read counts",
                count_of(distinct, "distinct row", "distinct rows")
            )
        },
    );
    let f = if fit.concentration.is_some() {
        "F fitted"
    } else if f_measurable {
        "F absent: its estimate stopped within 1% of an end of its search"
    } else {
        "F absent: no locus of the motif length has two samples with reads"
    };
    format!(
        "motif length {period} at {reference_repeats} repeats: {loci} with reads, {pairs} with \
         reads {distinct}; {f}",
        loci = count_of(evidence.tracts_with_reads(), "locus", "loci"),
        pairs = count_of(
            evidence.rows_with_reads(),
            "(locus, sample) pair",
            "(locus, sample) pairs"
        ),
    )
}

/// **How many standard errors were withheld because too few loci stand behind them**, over every
/// fitted class: every number fitted for a class with fewer than [`LOCI_FOR_AN_ERROR`] loci with
/// reads, and, in the other classes, each allele-length share expected at fewer than
/// [`LOCI_FOR_A_SHARE_ERROR`] loci. A kind with nothing withheld is left out, and with nothing
/// withheld at all the line says so. `None` when no class carries standard errors (none was
/// fitted, or none was asked for).
pub(super) fn withheld_errors_line(outcomes: &[StratumOutcome]) -> Option<String> {
    let mut any_errors = false;
    let (mut in_thin_classes, mut thin_classes, mut on_rare_lengths) = (0_usize, 0_usize, 0_usize);
    for outcome in outcomes {
        let StratumOutcome::Fitted(fit) = outcome else {
            continue;
        };
        let Some(errors) = &fit.standard_errors else {
            continue;
        };
        any_errors = true;
        let slippage = errors
            .slippage
            .iter()
            .flatten()
            .flat_map(|group| [group.level, group.shorter_share, group.fall_off]);
        let withheld = slippage
            .chain(errors.allele_length_variation.iter().copied())
            .chain([errors.concentration])
            .filter(|error| *error == StratumError::TooFewLoci)
            .count();
        if fit.tracts_of_its_own < LOCI_FOR_AN_ERROR {
            in_thin_classes += withheld;
            thin_classes += 1;
        } else {
            on_rare_lengths += withheld;
        }
    }
    any_errors.then(|| {
        let total = in_thin_classes + on_rare_lengths;
        if total == 0 {
            return "no standard error was withheld for too few loci".to_owned();
        }
        let mut parts = Vec::new();
        if thin_classes > 0 {
            parts.push(format!(
                "{in_thin_classes} in the {} with fewer than {LOCI_FOR_AN_ERROR} loci with reads \
                 (every number fitted for such a class)",
                count_of(thin_classes, "class", "classes"),
            ));
        }
        if on_rare_lengths > 0 {
            parts.push(format!(
                "{on_rare_lengths} on the shares of allele lengths expected at fewer than \
                 {LOCI_FOR_A_SHARE_ERROR} loci"
            ));
        }
        format!(
            "standard errors withheld because too few loci stand behind them: {total} in all; {}",
            parts.join(", ")
        )
    })
}

#[cfg(test)]
mod tests {
    use super::super::tests::{a_drawn_class, a_narrow_production_fit, a_stiff_three_base_cohort};
    use super::*;
    use crate::parameter_estimation::joint::ssr_fit::bench_fixtures::{draw_stratum, spectrum_of};
    use crate::parameter_estimation::joint::ssr_fit::expected_counts::DistinctRows;

    /// The fit's outcomes and every line it logged.
    fn fit_and_log(
        strata: &[StratumEvidence],
        homozygote_excess: &[f64],
    ) -> (Vec<StratumOutcome>, Vec<Line>) {
        let (outcomes, lines, _) = fit_log_and_rounds(strata, homozygote_excess);
        (outcomes, lines)
    }

    /// The fit's outcomes, every line it logged, and each round it reported as ended: its motif
    /// length and the rounds run so far.
    fn fit_log_and_rounds(
        strata: &[StratumEvidence],
        homozygote_excess: &[f64],
    ) -> (Vec<StratumOutcome>, Vec<Line>, Vec<(u8, usize)>) {
        let mut lines = Vec::new();
        let rounds = std::cell::RefCell::new(Vec::new());
        let outcomes = fit_repeat_tracts_logging(
            strata,
            homozygote_excess,
            &a_narrow_production_fit(),
            &mut |line| lines.push(line),
            &|period, round, _| rounds.borrow_mut().push((period, round)),
        );
        (outcomes, lines, rounds.into_inner())
    }

    /// The one line that starts with `prefix`.
    fn the_line<'a>(lines: &'a [Line], prefix: &str) -> &'a Line {
        let matching: Vec<&Line> = lines
            .iter()
            .filter(|line| line.text().starts_with(prefix))
            .collect();
        assert_eq!(
            matching.len(),
            1,
            "one line starting {prefix:?} in {lines:#?}"
        );
        matching[0]
    }

    /// **The log of a small drawn cohort holds every item of spec §2.9**, and logging changes no
    /// number.
    ///
    /// Two-base runs at 8, 9 and 10 repeats with 40, 5 and 40 loci, and one at 11 repeats with no
    /// reads; one three-base class with no reads, the only one of its motif length. Twelve samples
    /// at 8 reads. The class at 9 repeats is too thin for any error, so some are withheld.
    #[test]
    fn the_log_holds_every_item_the_spec_names() {
        let strata = vec![
            a_drawn_class(2, 8, 40, 12, 8, 31),
            a_drawn_class(2, 9, 5, 12, 8, 32),
            a_drawn_class(2, 10, 40, 12, 8, 33),
            a_drawn_class(2, 11, 0, 12, 8, 34),
            a_drawn_class(3, 7, 0, 12, 8, 35),
        ];
        let excess = [0.3; 12];
        let (outcomes, lines, rounds) = fit_log_and_rounds(&strata, &excess);
        assert_eq!(
            outcomes,
            fit_repeat_tracts(&strata, &excess, &a_narrow_production_fit()),
            "logging changes no number"
        );

        // How the fit ran: rounds, how it stopped, F, the time in each step. It is the one line
        // printed with time and memory.
        let ran = the_line(&lines, "motif length 2: ");
        assert!(matches!(ran, Line::StepEnded(_)), "{ran:?}");
        let ran = ran.text();
        for item in [
            "3 of 4 classes fitted together in ",
            " rounds, ",
            "F (the between-locus share: how much the allele-length frequencies vary from one \
             locus to the next) read off pairs of samples at a locus ",
            "time weighing genotypes and counting ",
            "drawing the curves ",
            "that also chose how stiff they are",
            "reading F off pairs ",
            "standard errors ",
        ] {
            assert!(ran.contains(item), "{item:?} missing from {ran}");
        }
        assert!(ran.contains("settled") || ran.contains("used all"), "{ran}");

        // Every round was reported as it ended, numbered on from the first, the last being the
        // rounds the line gives; the motif length with no reads ran none.
        let total: usize = ran
            .split(" fitted together in ")
            .nth(1)
            .and_then(|rest| rest.split(' ').next())
            .and_then(|count| count.parse().ok())
            .unwrap_or_else(|| panic!("the rounds in {ran}"));
        let expected: Vec<(u8, usize)> = (1..=total).map(|round| (2, round)).collect();
        assert_eq!(rounds, expected);

        // The five strengths, against the range searched.
        let stiffness = the_line(&lines, "motif length 2, how stiff each curve came out").text();
        for item in [
            "(searched from 1e-4, where each repeat count stands on its own, to 1e10, a straight \
             line): ",
            "slip rate ",
            "shorter share ",
            "fall-off ",
            "F ",
            "allele-length shares (one strength for all of them) ",
        ] {
            assert!(
                stiffness.contains(item),
                "{item:?} missing from {stiffness}"
            );
        }
        assert!(!stiffness.contains("F absent"), "{stiffness}");

        // Where bends were eased, or that none was.
        the_line(&lines, "motif length 2, bends eased");

        // One line a class: pairs with reads, distinct rows, F.
        for evidence in &strata[..3] {
            let repeats = evidence.stratum.reference_repeats;
            let line = the_line(&lines, &format!("motif length 2 at {repeats} repeats: ")).text();
            let distinct = DistinctRows::of(evidence, &excess).distinct_rows();
            let expected = format!(
                "{} loci with reads, {} (locus, sample) pairs with reads holding {distinct} \
                 distinct rows of read counts; F ",
                evidence.tracts_with_reads(),
                evidence.rows_with_reads()
            );
            assert!(line.contains(&expected), "{expected:?} missing from {line}");
            assert!(distinct < evidence.rows_with_reads(), "{line}");
        }
        assert_eq!(
            the_line(&lines, "motif length 2 at 11 repeats: ").text(),
            "motif length 2 at 11 repeats: no reads at any locus, so not fitted"
        );
        assert_eq!(
            the_line(&lines, "motif length 3: ").text(),
            "motif length 3: no reads at any locus in its 1 class, so nothing was fitted"
        );

        // The withheld errors, counted from the outcomes.
        let thin = match &outcomes[1] {
            StratumOutcome::Fitted(fit) => fit.as_ref(),
            other => panic!("the thin class has reads: {other:?}"),
        };
        let thin_withheld =
            3 + thin.allele_length_variation.len() + usize::from(thin.concentration.is_some());
        let withheld = lines
            .last()
            .map(Line::text)
            .filter(|line| line.contains("standard errors withheld"))
            .unwrap_or_else(|| panic!("the last line counts the withheld errors: {lines:#?}"));
        assert!(
            withheld.contains(&format!(
                "; {thin_withheld} in the 1 class with fewer than 8 loci with reads (every number \
                 fitted for such a class), "
            )),
            "{withheld}"
        );
        let on_rare_lengths: usize = [&outcomes[0], &outcomes[2]]
            .into_iter()
            .map(|outcome| match outcome {
                StratumOutcome::Fitted(fit) => fit
                    .standard_errors
                    .as_ref()
                    .expect("production asks for errors")
                    .allele_length_variation
                    .iter()
                    .filter(|error| **error == StratumError::TooFewLoci)
                    .count(),
                other => panic!("a rich class has reads: {other:?}"),
            })
            .sum();
        assert!(on_rare_lengths > 0, "the window's outer lengths are rare");
        assert!(
            withheld.starts_with(&format!(
                "standard errors withheld because too few loci stand behind them: {} in all; ",
                thin_withheld + on_rare_lengths
            )),
            "{withheld}"
        );
        assert!(
            withheld.ends_with(&format!(
                "{on_rare_lengths} on the shares of allele lengths expected at fewer than 3 loci"
            )),
            "{withheld}"
        );
    }

    /// **At one sample every class says F is absent, and why**: no locus holds two samples.
    #[test]
    fn at_one_sample_every_class_says_why_f_is_absent() {
        let strata: Vec<StratumEvidence> = (8..=10)
            .map(|repeats| a_drawn_class(1, repeats, 40, 1, 10, 20 + repeats))
            .collect();
        let (_, lines) = fit_and_log(&strata, &[0.3]);
        let ran = the_line(&lines, "motif length 1: ").text();
        assert!(
            ran.contains(
                "F (the between-locus share: how much the allele-length frequencies vary from one \
                 locus to the next) absent: no locus has two samples with reads"
            ),
            "{ran}"
        );
        let stiffness = the_line(&lines, "motif length 1, how stiff").text();
        assert!(stiffness.contains("F absent"), "{stiffness}");
        for repeats in 8..=10 {
            let line = the_line(&lines, &format!("motif length 1 at {repeats} repeats: ")).text();
            assert!(
                line.ends_with("F absent: no locus of the motif length has two samples with reads"),
                "{line}"
            );
        }
    }

    /// **A class whose F stopped at an end of its search says so**: 40 samples at 20 reads over 40
    /// loci drawn at a concentration of 10⁶, as the conversion's own test of that case.
    #[test]
    fn a_class_whose_f_stopped_at_an_end_of_its_search_says_so() {
        let strata: Vec<StratumEvidence> = (8..=9)
            .map(|repeats| {
                let truth = Slippage {
                    level: 0.03,
                    shorter_share: 0.8,
                    fall_off: 0.25,
                };
                let mut evidence = draw_stratum(
                    truth,
                    &spectrum_of(3),
                    1e6,
                    0.0,
                    40,
                    40,
                    20,
                    1,
                    50 + repeats,
                );
                evidence.stratum = Stratum {
                    period: 1,
                    reference_repeats: repeats,
                };
                evidence
            })
            .collect();
        let (_, lines) = fit_and_log(&strata, &[0.0; 40]);
        for repeats in 8..=9 {
            let line = the_line(&lines, &format!("motif length 1 at {repeats} repeats: ")).text();
            assert!(
                line.ends_with("F absent: its estimate stopped within 1% of an end of its search"),
                "{line}"
            );
        }
    }

    /// **When nothing is handed in, one line says there was nothing to fit**, and nothing else is
    /// logged, not even the count of withheld errors.
    #[test]
    fn one_line_says_so_when_there_are_no_classes() {
        let (outcomes, lines) = fit_and_log(&[], &[]);
        assert!(outcomes.is_empty());
        assert_eq!(
            lines,
            vec![Line::Detail(
                "no repeat classes were handed in, so there was nothing to fit".to_owned()
            )]
        );
    }

    /// **With one or two repeat counts with reads, no strength is printed**: the curve through
    /// them is a line whatever its strength, so the line says there is no bend to penalise, and
    /// no bend is reported eased. Spec §4.1's one-class and two-class cohorts.
    #[test]
    fn one_or_two_repeat_counts_have_no_bend_to_penalise() {
        let one = vec![a_drawn_class(2, 8, 40, 12, 8, 71)];
        let two = vec![
            a_drawn_class(2, 8, 40, 12, 8, 72),
            a_drawn_class(2, 9, 40, 12, 8, 73),
            a_drawn_class(2, 10, 0, 12, 8, 74),
        ];
        for (strata, said, which) in [
            (
                one,
                "one repeat count with reads",
                "its one class fitted in ",
            ),
            (
                two,
                "only two repeat counts with reads",
                "2 of 3 classes fitted together in ",
            ),
        ] {
            let (_, lines) = fit_and_log(&strata, &[0.3; 12]);
            let ran = the_line(&lines, "motif length 2: ").text();
            assert!(ran.contains(which), "{which:?} missing from {ran}");
            // No class is thin, so the withheld line prints no count of zero: either the shares
            // withheld alone, or that none was.
            let withheld: Vec<&str> = lines
                .iter()
                .map(Line::text)
                .filter(|line| line.contains("withheld"))
                .collect();
            assert_eq!(withheld.len(), 1, "{lines:#?}");
            assert!(!withheld[0].contains(" 0 "), "{}", withheld[0]);
            assert_eq!(
                the_line(&lines, "motif length 2, how stiff").text(),
                format!(
                    "motif length 2, how stiff each curve came out: {said}, so no curve has a \
                     bend to penalise; each class keeps its own counts"
                )
            );
            assert!(
                lines
                    .iter()
                    .all(|line| !line.text().contains("bends eased")),
                "{lines:#?}"
            );
            assert!(
                lines.iter().all(|line| !line.text().contains("searched")),
                "{lines:#?}"
            );
        }
    }

    /// **A where-it-has-got-to line names the motif length, the round, the largest change and the
    /// point at which the fit settles.**
    #[test]
    fn the_round_line_says_where_the_fit_has_got_to() {
        assert_eq!(
            round_line(2, 140, 3.14e-5, &a_narrow_production_fit()),
            "motif length 2, round 140 of at most 2000: the largest change in any number was \
             3.1e-5 (it settles once every change is below 1e-7)"
        );
    }

    /// **Plural nouns follow the count**, and a strength is written with two significant figures,
    /// plainly from 0.01 to below 1,000 and as a power of ten outside.
    #[test]
    fn counts_and_strengths_read_as_words() {
        assert_eq!(count_of(1, "class", "classes"), "1 class");
        assert_eq!(count_of(0, "class", "classes"), "0 classes");
        assert_eq!(count_of(5, "class", "classes"), "5 classes");
        for (strength, written) in [
            (2.39, "2.4"),
            (13.6, "14"),
            (136.0, "140"),
            (0.0586, "0.059"),
            (0.00586, "5.9e-3"),
            (1.56e9, "1.6e9"),
            (1e10, "1.0e10"),
        ] {
            assert_eq!(readable(strength), written, "{strength}");
        }
    }

    /// **Eased bends are named by curve and repeat count**, an allele-length share by how many
    /// repeats its alleles lie from the reference length; with none eased the line says so.
    #[test]
    fn eased_bends_are_named_by_curve_and_repeat_count() {
        let curve = |eased_at: Vec<u64>| ParameterCurve {
            first_repeats: 8,
            values: vec![0.1, 0.2, 0.3, 0.4],
            smoothing: 10.0,
            eased_at,
            strengths_rejected: 0,
            searches_without_a_strength: 0,
            unfactorised_fits_kept: 0,
        };
        let mut curves = SmoothCurves {
            level: curve(vec![9]),
            shorter_share: curve(Vec::new()),
            fall_off: curve(Vec::new()),
            between_locus_share: None,
            allele_length_variation: vec![
                curve(Vec::new()),
                curve(vec![9, 10]),
                curve(Vec::new()),
                curve(vec![11]),
                curve(vec![10]),
            ],
        };
        assert_eq!(
            eased_bends_line(2, &curves, 2),
            "motif length 2, bends eased because the reads insisted on a kink: slip rate at 9 \
             repeats; share of alleles 1 repeat shorter at 9, 10 repeats; share of alleles 1 \
             repeat longer at 11 repeats; share of alleles 2 repeats longer at 10 repeats"
        );
        curves.level.eased_at.clear();
        for share in &mut curves.allele_length_variation {
            share.eased_at.clear();
        }
        assert_eq!(
            eased_bends_line(2, &curves, 2),
            "motif length 2, bends eased: none; no curve's reads insisted on a kink"
        );
    }

    /// **No line about withheld errors when no class carries errors**: a fit not asked for them.
    #[test]
    fn no_withheld_line_without_errors() {
        let strata = vec![a_drawn_class(1, 8, 10, 3, 5, 61)];
        let config = PenalisedFitConfig {
            standard_errors: false,
            ..a_narrow_production_fit()
        };
        let outcomes = fit_repeat_tracts(&strata, &[0.3; 3], &config);
        assert_eq!(withheld_errors_line(&outcomes), None);
    }

    /// **A motif length whose errors could not be determined says so, and writes every error as
    /// not placed** (or withheld for too few loci, which comes first); the stiff three-base cohort
    /// of the fit's own test, 10 loci a class.
    #[test]
    fn undetermined_errors_are_logged_and_written_as_not_placed() {
        let (outcomes, lines) = fit_and_log(&a_stiff_three_base_cohort(), &[0.3; 4]);
        let line = the_line(&lines, "motif length 3, standard errors:");
        assert_eq!(
            line.text(),
            "motif length 3, standard errors: none determined for any of its numbers except F, \
             which is read off pairs of samples with a curvature of its own: the fit's \
             curvature, with the curves' bend penalties added, could not be inverted, either \
             because some combination of numbers has next to no reads behind it (the slope of a \
             stiff curve through a single repeat count's slipped reads, for one) or because the \
             curvature is not positive in some direction at the values the fit reached; those \
             errors are written as not placed, except where too few loci already withheld them"
        );
        let mut errors = 0;
        for outcome in &outcomes {
            let StratumOutcome::Fitted(fit) = outcome else {
                panic!("every class has reads");
            };
            let every = fit.standard_errors.as_ref().expect("errors were asked for");
            let slippage = every
                .slippage
                .iter()
                .flatten()
                .flat_map(|group| [group.level, group.shorter_share, group.fall_off]);
            for error in slippage.chain(every.allele_length_variation.iter().copied()) {
                assert!(
                    matches!(error, StratumError::NotPlaced | StratumError::TooFewLoci),
                    "{error:?}"
                );
                errors += 1;
            }
            // F's own curvature, read off pairs of samples, was inverted: its error is kept.
            assert!(
                matches!(every.concentration, StratumError::Estimated(_)),
                "{:?}",
                every.concentration
            );
        }
        assert!(errors > 0);
    }

    /// **Strengths that could not be scored are named by curve**, and no line is printed when
    /// every strength tried was scored.
    #[test]
    fn unscored_strengths_are_named_by_curve() {
        let curve = |rejected: usize, kept: usize| ParameterCurve {
            first_repeats: 8,
            values: vec![0.1, 0.2, 0.3],
            smoothing: 10.0,
            eased_at: Vec::new(),
            strengths_rejected: rejected,
            searches_without_a_strength: kept,
            unfactorised_fits_kept: 0,
        };
        let mut curves = SmoothCurves {
            level: curve(0, 0),
            shorter_share: curve(0, 0),
            fall_off: curve(0, 0),
            between_locus_share: Some(curve(0, 0)),
            allele_length_variation: vec![curve(0, 0); 3],
        };
        assert_eq!(unscored_strengths_line(2, &curves), None);
        curves.shorter_share = curve(3, 0);
        curves.between_locus_share = Some(curve(5, 1));
        assert_eq!(
            unscored_strengths_line(2, &curves).as_deref(),
            Some(
                "motif length 2, strengths that could not be scored: strengths never chosen \
                 because the curve's fit there could not be factorised, even with a larger \
                 ridge: shorter share 3, F 5; no strength tried could be scored, so the strength \
                 held before was kept: F in 1 search"
            )
        );
        // A fit kept at a strength no search chose, alone, is enough for the line.
        let mut curves = SmoothCurves {
            level: curve(0, 0),
            shorter_share: curve(0, 0),
            fall_off: curve(0, 0),
            between_locus_share: None,
            allele_length_variation: vec![curve(0, 0); 3],
        };
        curves.fall_off.unfactorised_fits_kept = 2;
        assert_eq!(
            unscored_strengths_line(2, &curves).as_deref(),
            Some(
                "motif length 2, strengths that could not be scored: fits at a strength no search \
                 chose that could not be factorised, so the curve kept the values Newton's method \
                 had reached when it stopped: fall-off 2"
            )
        );
    }

    /// **The two variants that keep some errors name what they keep**: F's alone undetermined,
    /// and every number's but F's.
    #[test]
    fn undetermined_errors_lines_name_what_is_kept() {
        assert_eq!(
            undetermined_errors_line(2, UndeterminedErrors::BetweenLocusShare),
            "motif length 2, standard errors: none determined for F: the curvature of F read off \
             pairs of samples, with its curve's bend penalty added, could not be inverted, either \
             because some combination of numbers has next to no reads behind it (the slope of a \
             stiff curve through a single repeat count's slipped reads, for one) or because the \
             curvature is not positive in some direction at the values the fit reached; those \
             errors are written as not placed, except where too few loci already withheld them"
        );
        assert_eq!(
            undetermined_errors_line(2, UndeterminedErrors::EveryNumberButBetweenLocusShare),
            "motif length 2, standard errors: none determined for any of its numbers except F, \
             which is read off pairs of samples with a curvature of its own: the fit's curvature, \
             with the curves' bend penalties added, could not be inverted, either because some \
             combination of numbers has next to no reads behind it (the slope of a stiff curve \
             through a single repeat count's slipped reads, for one) or because the curvature is \
             not positive in some direction at the values the fit reached; those errors are \
             written as not placed, except where too few loci already withheld them"
        );
    }
}
