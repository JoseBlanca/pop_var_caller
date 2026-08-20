//! How the slippage level rises with repeat count, and where an emitted level came from.
//!
//! **The problem this exists for.** The four slippage numbers are fitted separately in every
//! (motif period, repeat count) cell, so neighbouring cells are free to disagree, and a cell too
//! thin to fit takes a neighbour's value whole — which costs 15 to 25% of the level per repeat
//! count borrowed across. On the tomato bench run 65 of 71 cells got no answer from their own
//! tracts at all.
//!
//! **What this module holds.** The *level* — how often a read reports a tract length other than
//! its allele's — becomes a curve in repeat count, fitted once per motif period and evaluated at
//! every cell. The other three numbers are untouched and keep their per-cell fit
//! (`doc/devel/ng/spec/str_slippage_level_curve.md` §1.2).
//!
//! **The family, and why it has a shape number instead of a fixed shape.** Two cohorts fitted
//! with nothing linking their cells prefer *opposite* shapes over the same repeat counts: over
//! 8 to 12 repeats, tomato's homopolymers are predicted best by an exponential (12.4% held-out
//! against a straight line's 33.6%) and HG002's by a straight line (8.0% against 31.2%). And both
//! fixed shapes produce impossible numbers a few repeat counts outside the range they saw — an
//! exponential fitted on HG002's 8-to-12 cells says the level at 30 repeats is 205, where the cell
//! at 30 repeats fits 0.120. So the shape is fitted rather than chosen
//! (`doc/devel/ng/reports/str_slippage_shape_2026-08-20.md` §4).
//!
//! Design: `doc/devel/ng/spec/str_slippage_level_curve.md`. Build order:
//! `doc/devel/ng/impl_plan/str_slippage_level_curve.md`.

use std::fmt;

// ---------------------------------------------------------------------
// The shape number
// ---------------------------------------------------------------------

/// How the slippage level compounds as repeat count rises.
///
/// **At 0 each extra repeat multiplies the level by a fixed factor; at 1 each extra repeat adds
/// a fixed amount.** In between it is the power the level is raised to before a straight line is
/// fitted through it:
///
/// ```text
/// level ^ rise_shape  =  intercept + slope · repeat_count      for rise_shape > 0
/// log(level)          =  intercept + slope · repeat_count      for rise_shape = 0
/// ```
///
/// **1 is production's own shape** — `baseline + slope · units`, clamped, in
/// `src/ssr/cohort/em.rs`'s `candidate_level` — so the family contains what production does
/// rather than replacing it with something unrelated.
///
/// Fitted per motif period and shared by every slippage group at that period, because a
/// curvature needs the whole span of repeat counts visible at once where a level and a slope do
/// not (spec §3).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct RiseShape(f64);

impl RiseShape {
    /// The multiplying end: each extra repeat multiplies the level by a fixed factor.
    pub const MULTIPLYING: Self = Self(0.0);

    /// The adding end: each extra repeat adds a fixed amount, which is production's shape.
    pub const ADDING: Self = Self(1.0);

    /// `None` outside `[0, 1]` or on a value that is not a number.
    ///
    /// **The range is closed at both ends and not widened**: past 1 the curve is convex in
    /// repeat count, which would say slippage accelerates without limit, and no cell of either
    /// cohort asks for it — HG002's homopolymers land exactly at 1 and its dinucleotides at 0.8.
    pub fn new(value: f64) -> Option<Self> {
        (value.is_finite() && (0.0..=1.0).contains(&value)).then_some(Self(value))
    }

    pub fn get(self) -> f64 {
        self.0
    }

    /// Whether this is the multiplying end, where the straight line is fitted to the logarithm.
    ///
    /// **The two branches are genuinely different arithmetic**, not a limit taken numerically:
    /// the power transform `(level^s − 1) / s` is undefined at `s = 0`, and its limit there is
    /// `log(level)`.
    pub fn is_multiplying(self) -> bool {
        self.0 == 0.0
    }
}

impl fmt::Display for RiseShape {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:.2}", self.0)
    }
}

// ---------------------------------------------------------------------
// What the curve is fitted from, and what it is
// ---------------------------------------------------------------------

/// One cell's own answer, as the curve fit sees it.
///
/// **Only cells that were fitted from their own tracts are here.** A cell that borrowed, or that
/// was refused, has nothing to say about the shape and would pull hardest if it were let in
/// (spec §4.1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FittedCell {
    /// The reference tract's repeat count — the axis the curve runs along.
    pub repeats: u64,
    /// The cell's own fitted slippage level, in `(0, 1)`.
    pub level: f64,
    /// **How many of the cell's reads the fitted level says slipped** — `level × reads crossing`.
    /// This is what sets how precisely the cell determines its own level, and it is the weight
    /// the curve fit uses (spec §4.2).
    ///
    /// **Not the count of reads sitting off the reference tract length**, which is much larger
    /// and at a polymorphic tract is mostly genuine allele length: at HG002's 30-repeat
    /// homopolymer cell, 60 reads in 100 that cross the tract report a length other than the
    /// reference's and the fit attributes 12 of those 60 to slippage.
    pub slipped_reads: f64,
}

impl FittedCell {
    /// How precisely this cell determines its own level, as a share of the level itself.
    ///
    /// `1 / sqrt(slipped reads)` — the best case, since the level is a proportion over the reads
    /// that moved. The real estimator sums over the genotype and is noisier than that, so this
    /// **understates** the cell's error and therefore leans the blend of spec §7 toward the cell.
    ///
    /// **Fewer than one slipped read is read as one**, so the answer is a share of the level and
    /// never infinite. Such a cell determines nothing and does not feed a curve (spec §4.1); the
    /// clamp exists so that weighing it against a curve is arithmetic rather than a special case.
    pub fn relative_standard_error(&self) -> f64 {
        1.0 / self.slipped_reads.max(1.0).sqrt()
    }
}

/// How the slippage level rises with repeat count, for one slippage group at one motif period.
///
/// **Held flat outside the repeat counts it was fitted over, never continued** — see
/// [`SlippageCurve::level_at`] and spec §6.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlippageCurve {
    /// How the level compounds; shared by every slippage group at this motif period.
    pub rise_shape: RiseShape,
    /// The straight line fitted to `level ^ rise_shape` against repeat count.
    pub intercept: f64,
    /// The same line's slope. **Positive by construction**: a fit returning a level that falls
    /// with repeat count is refused rather than emitted (spec §9).
    pub slope: f64,
    /// The lowest repeat count of a cell that fed this curve.
    pub fitted_from: u64,
    /// The highest. Beyond `fitted_from ..= fitted_to` the level is held at the nearer end.
    pub fitted_to: u64,
    /// **How far the curve landed from a cell it had not seen** — the median relative error over
    /// leaving each contributing cell out in turn. Spec §7 reads it as the curve's own precision
    /// when weighing the curve against a cell's own answer.
    pub held_out_error: f64,
    /// How many cells stood behind it. A curve through four cells and one through twenty-three
    /// are both curves, and a consumer must be able to tell them apart.
    pub cells: usize,
}

impl SlippageCurve {
    /// The slippage level this curve gives at `repeats`.
    ///
    /// **Inside the repeat counts the curve was fitted over it is the fitted line, read back
    /// through the shape number. Outside, it is held at the nearer fitted end.**
    ///
    /// Holding flat rather than continuing the line is the whole of spec §6, and the reason is
    /// that continuing it produces numbers that are not probabilities: an exponential fitted on
    /// HG002's 8-to-12 homopolymer cells says 205 at 30 repeats, where the cell at 30 repeats
    /// fits 0.120. **Held flat is wrong in a known direction** — slippage genuinely keeps rising,
    /// so a tract longer than anything fitted is under-stated — and [`SlippageCurve::reach`] is
    /// what lets a consumer see it.
    pub fn level_at(&self, repeats: u64) -> f64 {
        assert!(
            self.fitted_from <= self.fitted_to,
            "a curve fitted from {} to {} repeats has no range to hold at its ends",
            self.fitted_from,
            self.fitted_to
        );
        let inside = repeats.clamp(self.fitted_from, self.fitted_to);
        self.level_on_the_line(inside as f64)
    }

    /// Whether `repeats` sat inside the repeat counts this curve was fitted over.
    pub fn reach(&self, repeats: u64) -> CurveReach {
        if repeats < self.fitted_from {
            CurveReach::BelowFitted
        } else if repeats > self.fitted_to {
            CurveReach::AboveFitted
        } else {
            CurveReach::Inside
        }
    }

    /// The fitted line read back into a level, with no regard for the fitted range.
    ///
    /// **Clamped into `(0, 1)` at both ends**, because the line is fitted in a transformed space
    /// where nothing stops it crossing either boundary: at the multiplying end the exponential
    /// passes 1 a few repeat counts above its range, and at the adding end the line passes 0 a
    /// few below — HG002's homopolymer line goes negative below 7.4 repeats. The clamp is not
    /// what makes extrapolation safe; holding at the fitted ends is. It is here so that no
    /// caller can be handed a number that is not a probability.
    fn level_on_the_line(&self, repeats: f64) -> f64 {
        let line = self.intercept + self.slope * repeats;
        let level = if self.rise_shape.is_multiplying() {
            line.exp()
        } else {
            // `level ^ s = line`, so `level = line ^ (1/s)`; a line at or below zero has no
            // real root and means the level has fallen off the bottom of the family.
            if line <= 0.0 {
                0.0
            } else {
                line.powf(1.0 / self.rise_shape.get())
            }
        };
        if level.is_finite() {
            level.clamp(LEVEL_FLOOR, LEVEL_CEILING)
        } else {
            LEVEL_CEILING
        }
    }
}

/// The smallest level a curve may report.
///
/// **Not zero.** A level of exactly zero says a read can never misread the tract, which the
/// emission model then treats as certainty; and the blend of spec §7 works on the logarithm,
/// which has no value there. One read in ten million is far below anything either cohort
/// measures — the thinnest cell fitted is 3.7 in 1,000 — so the floor cannot bind on real data
/// and exists to keep the arithmetic total.
pub const LEVEL_FLOOR: f64 = 1e-7;

/// The largest level a curve may report.
///
/// **Not one.** At a level of one every read has misread, and the emission model divides by the
/// share of reads that did not. HG002's longest homopolymer cell fits 0.120, so a curve reaching
/// this has been asked for a repeat count far outside anything fitted.
pub const LEVEL_CEILING: f64 = 0.999;

// ---------------------------------------------------------------------
// Fitting one line, at a shape number someone else chose
// ---------------------------------------------------------------------

/// Why a set of cells produced no curve.
///
/// **No `Eq`**: one variant carries the slope it refused, which is a float.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NoCurve {
    /// Fewer usable cells than a line needs.
    ///
    /// **Two floors exist and they are different questions.** [`fit_line`] refuses below **two**,
    /// which is arithmetic — one point does not make a line. [`choose_rise_shape`] refuses below
    /// [`SlippageCurveConfig::min_cells_for_a_curve`], which is about whether the answer means
    /// anything. `floor` says which one refused.
    TooFewCells { cells: usize, floor: usize },
    /// Every cell sits at the same repeat count, so no line through them has a slope.
    OneRepeatCountOnly { repeats: u64 },
    /// The best line through these cells says slippage *falls* as tracts get longer.
    ///
    /// **Refused rather than emitted.** Every measurement this design rests on says the level
    /// rises with repeat count, so a falling fit is reporting the cells' noise; a consumer handed
    /// it would give a long tract a smaller error rate than a short one (spec §9).
    LevelWouldFall { slope: f64 },
}

/// Fit one slippage group's line through its cells, at a shape number chosen elsewhere.
///
/// **Weighted least squares of `level ^ rise_shape` on repeat count**, each cell weighted by its
/// slipped reads — the count that sets how precisely that cell determines its own level. The
/// `held_out_error` and `cells` of the returned curve are placeholders that
/// [`choose_rise_shape`] fills; a caller that fits one line directly gets the cell count and a
/// held-out error of zero, and must not read the second as a claim.
///
/// **Cells whose level is not strictly inside `(0, 1)` are dropped before fitting.** The level is
/// bounded below by zero and a thin cell piles up against that boundary, where the multiplying
/// end's logarithm has no value (spec §7's last note).
pub fn fit_line(cells: &[FittedCell], rise_shape: RiseShape) -> Result<SlippageCurve, NoCurve> {
    let usable: Vec<&FittedCell> = cells
        .iter()
        .filter(|cell| cell.level > 0.0 && cell.level < 1.0 && cell.slipped_reads > 0.0)
        .collect();
    if usable.len() < 2 {
        return Err(NoCurve::TooFewCells {
            cells: usable.len(),
            floor: 2,
        });
    }

    let transform = |level: f64| {
        if rise_shape.is_multiplying() {
            level.ln()
        } else {
            level.powf(rise_shape.get())
        }
    };

    // Weighted means first, so the two sums below are about the centred values and cannot lose
    // precision to a large common offset.
    let total_weight: f64 = usable.iter().map(|cell| cell.slipped_reads).sum();
    let mean_repeats: f64 = usable
        .iter()
        .map(|cell| cell.slipped_reads * cell.repeats as f64)
        .sum::<f64>()
        / total_weight;
    let mean_transformed: f64 = usable
        .iter()
        .map(|cell| cell.slipped_reads * transform(cell.level))
        .sum::<f64>()
        / total_weight;

    let mut covariance = 0.0;
    let mut spread = 0.0;
    for cell in &usable {
        let from_mean = cell.repeats as f64 - mean_repeats;
        covariance += cell.slipped_reads * from_mean * (transform(cell.level) - mean_transformed);
        spread += cell.slipped_reads * from_mean * from_mean;
    }
    if spread <= 0.0 {
        return Err(NoCurve::OneRepeatCountOnly {
            repeats: usable[0].repeats,
        });
    }

    // **Not `slope <= 0.0`.** That comparison is false for a slope that is not a number, which
    // would let a `NaN` line through and hand every cell a `NaN` level.
    let slope = covariance / spread;
    if !(slope.is_finite() && slope > 0.0) {
        return Err(NoCurve::LevelWouldFall { slope });
    }

    Ok(SlippageCurve {
        rise_shape,
        intercept: mean_transformed - slope * mean_repeats,
        slope,
        fitted_from: usable
            .iter()
            .map(|cell| cell.repeats)
            .min()
            .expect("two cells"),
        fitted_to: usable
            .iter()
            .map(|cell| cell.repeats)
            .max()
            .expect("two cells"),
        held_out_error: 0.0,
        cells: usable.len(),
    })
}

// ---------------------------------------------------------------------
// Choosing the shape number, across the slippage groups of one period
// ---------------------------------------------------------------------

/// One motif period's curves: one line per slippage group, all sharing a shape number.
///
/// **The shape number is shared and the lines are not**, because they answer different questions
/// with different amounts of evidence behind them. A line's level and slope are where a library's
/// own chemistry lives and a period with four cells can carry them; the curvature needs the whole
/// span of repeat counts visible at once, and one library of a 63-accession cohort at three reads
/// a position puts about twelve slipped reads behind a cell (spec §3).
#[derive(Debug, Clone, PartialEq)]
pub struct PeriodCurves {
    /// The shape number every group at this period shares.
    pub rise_shape: RiseShape,
    /// One curve a slippage group, indexed by group, `None` where that group had no line.
    pub by_group: Vec<Option<SlippageCurve>>,
    /// How far the winning shape number landed from cells it had not seen — the median relative
    /// error over leaving each contributing cell of every group out in turn.
    pub held_out_error: f64,
    /// Contributing cells, over every group.
    pub cells: usize,
}

/// Choose the shape number for one motif period, and fit every slippage group's line at it.
///
/// `cells_by_group` is one list of contributing cells per slippage group, indexed the way
/// `StratumFit::slippage` is; a group with no cells contributes nothing and gets `None` back.
///
/// **The rung is chosen on cells it did not see.** For each rung of the grid, every group's line
/// is fitted, then each contributing cell is left out in turn, its group's line refitted without
/// it, and the cell predicted; the rung with the lowest median relative error wins. Fit quality
/// on the cells a rung *did* see decides nothing — any rung can be made to fit what it saw, and
/// the adding end has no more freedom than the multiplying one, so only held-out prediction
/// separates them (spec §4.3).
///
/// **Ties go to the larger rung.** The two ends fail differently: a level that is too small
/// under-states real slippage, where the multiplying end returns numbers above one a few repeat
/// counts outside its range and has to be clamped.
pub fn choose_rise_shape(
    cells_by_group: &[Vec<FittedCell>],
    config: &SlippageCurveConfig,
) -> Result<PeriodCurves, NoCurve> {
    let contributing: usize = cells_by_group.iter().map(Vec::len).sum();
    if contributing < config.min_cells_for_a_curve {
        return Err(NoCurve::TooFewCells {
            cells: contributing,
            floor: config.min_cells_for_a_curve,
        });
    }

    let mut best: Option<(f64, PeriodCurves)> = None;
    for rise_shape in config.rise_shape_grid() {
        let by_group: Vec<Option<SlippageCurve>> = cells_by_group
            .iter()
            .map(|cells| fit_line(cells, rise_shape).ok())
            .collect();
        if by_group.iter().all(Option::is_none) {
            continue;
        }
        let Some(held_out_error) = held_out_error_of(cells_by_group, rise_shape) else {
            continue;
        };
        let fitted: usize = by_group
            .iter()
            .flatten()
            .map(|curve| curve.cells)
            .sum::<usize>();
        let curves = PeriodCurves {
            rise_shape,
            by_group,
            held_out_error,
            cells: fitted,
        };
        // `<=` rather than `<`: the grid runs from the multiplying end upward, so an equal score
        // at a later rung replaces an earlier one, which is the tie rule.
        if best
            .as_ref()
            .is_none_or(|(best_error, _)| held_out_error <= *best_error)
        {
            best = Some((held_out_error, curves));
        }
    }

    let Some((held_out_error, mut curves)) = best else {
        return Err(NoCurve::LevelWouldFall { slope: f64::NAN });
    };
    // Every emitted curve carries the period's own held-out error and cell count, so a consumer
    // holding one curve can weigh it without holding the period.
    for curve in curves.by_group.iter_mut().flatten() {
        curve.held_out_error = held_out_error;
        curve.cells = curves.cells;
    }
    Ok(curves)
}

/// The median relative error of leaving each contributing cell out in turn, at one rung.
///
/// **The score continues the line where a deployed curve holds flat, and the difference is
/// deliberate.** Leaving out the lowest or highest cell puts it outside what remains, so a
/// deployed curve would predict it with its neighbour's value — which says nothing about the
/// shape, and the end cells are exactly where a shape shows. Measured three ways on both
/// cohorts, the chosen rung is the same at four of the five periods that can be scored; at
/// HG002's dinucleotides, holding flat instead moves it from 0.80 to 0.70, and scoring only the
/// cells that stay inside the range leaves 3 of tomato's 5 cells scored and moves its answer
/// from 0.00 to 0.15. Continuing the line is the variant that keeps every cell in the score.
///
/// **This does not leak into what a caller reads.** [`SlippageCurve::level_at`] still holds flat;
/// only the rung's score continues the line, and it is thrown away once the rung is picked.
///
/// `None` when no cell could be predicted at all — every leave-one-out fit refused, which
/// happens when a group has two cells and dropping one leaves a point.
fn held_out_error_of(cells_by_group: &[Vec<FittedCell>], rise_shape: RiseShape) -> Option<f64> {
    let mut errors: Vec<f64> = Vec::new();
    for cells in cells_by_group {
        for held in 0..cells.len() {
            let rest: Vec<FittedCell> = cells
                .iter()
                .enumerate()
                .filter(|(index, _)| *index != held)
                .map(|(_, cell)| *cell)
                .collect();
            let Ok(curve) = fit_line(&rest, rise_shape) else {
                continue;
            };
            let target = cells[held];
            if target.level <= 0.0 {
                continue;
            }
            let predicted = curve.level_on_the_line(target.repeats as f64);
            errors.push(((predicted - target.level) / target.level).abs());
        }
    }
    if errors.is_empty() {
        return None;
    }
    errors.sort_by(|left, right| left.total_cmp(right));
    let middle = errors.len() / 2;
    Some(if errors.len().is_multiple_of(2) {
        (errors[middle - 1] + errors[middle]) / 2.0
    } else {
        errors[middle]
    })
}

// ---------------------------------------------------------------------
// Where an emitted level came from
// ---------------------------------------------------------------------

/// Where a cell's emitted slippage level came from.
///
/// **After this change a level fitted from 8,000 slipped reads and one interpolated across a gap
/// look the same in the number alone**, and the mechanism that used to mark the second — the
/// borrowing rule — no longer sets the level. This is its replacement (spec §8).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LevelSource {
    /// The cell's own fit, whole. Its period had no curve.
    Cell,
    /// The curve, whole. The cell had no fit of its own.
    Curve,
    /// Both, weighted by how precisely each determines the level.
    ///
    /// `curve_weight` is the share the curve carried, in `[0, 1]`. It runs about 0.93 at a cell
    /// with 40 slipped reads behind it and about 0.06 at one with 8,000, at a curve whose
    /// held-out error is 4.4% (spec §7).
    Blend { curve_weight: f64 },
}

impl LevelSource {
    /// The share the curve carried, whichever variant this is — 1 for [`LevelSource::Curve`] and
    /// 0 for [`LevelSource::Cell`].
    ///
    /// **The three variants are one formula with the weight at its two ends and in between**, so
    /// a consumer that only wants the number should not have to match on which end it is.
    pub fn curve_weight(self) -> f64 {
        match self {
            Self::Cell => 0.0,
            Self::Curve => 1.0,
            Self::Blend { curve_weight } => curve_weight,
        }
    }
}

/// What one cell's emitted slippage level is, once its own fit and its period's curve have both
/// had their say.
///
/// **Inverse-variance on the log scale.** Both quantities are relative errors, so they combine
/// multiplicatively: the cell's own error is `1 / sqrt(slipped reads)` and the curve's is its
/// held-out error, and each gets weight in proportion to `1 / error²`.
///
/// **This is the protection against fitting each cell's noise, and it is not a switch between
/// two regimes.** At a curve whose held-out error is 4.4%, the curve carries about 93% of the
/// weight at a cell with 40 slipped reads and about 6% at one with 8,000. Using the curve
/// everywhere is this formula with the curve's weight pinned at one; using each cell's own
/// answer is it pinned at zero.
///
/// **Why not pin it at one.** Where the curve misses it misses systematically, and at the cells
/// holding the most tracts: over HG002's 23 homopolymer cells the winning curve sits within 0.5
/// to 12% of every cell from 10 repeats up and is 27% and 55% high at 8 and 9, against those
/// cells' own errors of 1.8% and 1.7%. Always-curve would report 10.4 reads slipping per 1,000
/// at a 9-repeat homopolymer where that cell's own 3,520 slipped reads say 6.7
/// (`str_slippage_level_curve.md` §7.1).
///
/// **The knee, and it is a refinement rather than a rescue.** A disagreement far larger than
/// either error explains is evidence about the curve, not about the cell, so beyond
/// `disagreement_knee` combined errors the curve's weight is divided by `(gap / knee)²`. At that
/// 9-repeat cell the gap is 9.3 combined errors, which no sampling noise produces. Without the
/// knee the blend already emits 7.1 there rather than 10.4, because a cell with 3,520 slipped
/// reads outweighs the curve on its own; the knee takes it to 6.8, and it is worth the one
/// comparison because the bottom of the repeat range is what the copy-floor decision reads.
pub fn blend_level(
    cell: Option<FittedCell>,
    curve: Option<&SlippageCurve>,
    repeats: u64,
    config: &SlippageCurveConfig,
) -> Option<BlendedLevel> {
    match (cell, curve) {
        // No curve at this period: the cell keeps its own answer, exactly as today.
        (Some(cell), None) => Some(BlendedLevel {
            level: cell.level,
            source: LevelSource::Cell,
            reach: None,
        }),
        // No fit of the cell's own: the curve supplies the level whole. This is the case the
        // whole change exists for — a stratum below the refusal floor used to get nothing.
        (None, Some(curve)) => Some(BlendedLevel {
            level: curve.level_at(repeats),
            source: LevelSource::Curve,
            reach: Some(curve.reach(repeats)),
        }),
        (Some(cell), Some(curve)) => {
            let from_curve = curve.level_at(repeats);
            // **A cell whose level fitted to zero has no logarithm and no evidence.** It takes
            // the curve whole rather than dragging the blend to zero.
            if cell.level <= 0.0 || cell.slipped_reads <= 0.0 {
                return Some(BlendedLevel {
                    level: from_curve,
                    source: LevelSource::Curve,
                    reach: Some(curve.reach(repeats)),
                });
            }
            let cell_error = cell.relative_standard_error();
            let curve_error = curve.held_out_error.max(f64::EPSILON);

            let gap = (cell.level.ln() - from_curve.ln()).abs()
                / (cell_error * cell_error + curve_error * curve_error).sqrt();
            let trust = if gap > config.disagreement_knee {
                let over = gap / config.disagreement_knee;
                1.0 / (over * over)
            } else {
                1.0
            };

            let cell_weight = 1.0 / (cell_error * cell_error);
            let curve_weight = trust / (curve_error * curve_error);
            let total = cell_weight + curve_weight;
            let share_of_curve = curve_weight / total;
            Some(BlendedLevel {
                level: ((1.0 - share_of_curve) * cell.level.ln()
                    + share_of_curve * from_curve.ln())
                .exp()
                .clamp(LEVEL_FLOOR, LEVEL_CEILING),
                source: LevelSource::Blend {
                    curve_weight: share_of_curve,
                },
                reach: Some(curve.reach(repeats)),
            })
        }
        // Neither: nothing to emit, and saying so is the difference between missing and quiet.
        (None, None) => None,
    }
}

/// One cell's emitted level and where it came from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BlendedLevel {
    pub level: f64,
    pub source: LevelSource,
    /// Whether the cell sat inside the curve's fitted range; `None` where there was no curve.
    pub reach: Option<CurveReach>,
}

/// Whether a cell sat inside the repeat counts its curve was fitted over.
///
/// **Held flat is wrong in a known direction and that is why it is recorded**: the level
/// genuinely keeps rising, so a cell above the fitted range is under-stated (spec §6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveReach {
    /// Inside `fitted_from ..= fitted_to`.
    Inside,
    /// Below the lowest cell that fed the curve; the level is that cell's.
    BelowFitted,
    /// Above the highest; the level is that cell's.
    AboveFitted,
}

// ---------------------------------------------------------------------
// The knobs, and where each number comes from
// ---------------------------------------------------------------------

/// How many rungs the shape number is searched over, spanning `[0, 1]` inclusive.
///
/// **21, which is steps of 0.05.** The two cohorts' periods land at 0.00, 0.80 and 1.00, so the
/// grid has to resolve 0.05 to distinguish period 2's answer from the adding end; nothing
/// measured asks for finer, and each extra rung costs one weighted least-squares fit per cell per
/// slippage group.
pub const RISE_SHAPE_RUNGS: usize = 21;

/// How few contributing cells leave a period without a curve.
///
/// **4, and it is arithmetic rather than a measurement — the smallest count at which leaving one
/// cell out still leaves a line and a spare.** It is soft, and spec §11 records the measurement
/// that would settle it along with the reason to expect it is too low: HG002's period 3 has
/// exactly four cells and its best rung predicts a held-out cell only to 31%, against 3.8% at
/// period 2's twenty.
pub const MIN_CELLS_FOR_A_CURVE: usize = 4;

/// How far a cell may sit from the curve, in the two errors combined, before the curve is taken
/// to be wrong about that cell rather than the cell unlucky.
///
/// **2.5 combined errors.** The case it exists for is measured: at HG002's 9-repeat homopolymer
/// the curve is 55% high against a cell whose own sampling error is 1.7%, a gap of **9.3**
/// combined errors, which no sampling noise produces. *The value is a conventional outlier knee
/// and is soft* — spec §7.2 records that without it the blend was already within 5.8% there
/// rather than 55%, because a cell with 3,520 slipped reads outweighs the curve on its own.
pub const DISAGREEMENT_KNEE: f64 = 2.5;

/// What a run may change about how the level's curve is fitted.
///
/// **The weight each cell carries is deliberately not here.** It is the cell's slipped reads, and
/// the measurement says the choice barely matters: the winning family's held-out error moves from
/// 5.13% unweighted to 4.39% weighted by slipped reads, and the ranking of families is identical
/// under all four weights tried. A knob nobody sets is how a measured decision gets changed by
/// accident.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SlippageCurveConfig {
    /// Rungs of the shape-number grid; see [`RISE_SHAPE_RUNGS`].
    pub rise_shape_rungs: usize,
    /// Below this many contributing cells a period gets no curve; see [`MIN_CELLS_FOR_A_CURVE`].
    pub min_cells_for_a_curve: usize,
    /// Where the curve stops being trusted about a cell; see [`DISAGREEMENT_KNEE`].
    pub disagreement_knee: f64,
    /// **`false` keeps every cell's own fitted level and draws no curve at all**, which is the
    /// arm the parity oracle runs: nothing here touches how a cell is fitted, so a cell's own
    /// level moving with the curve off is a defect in the plumbing.
    pub draw_curves: bool,
}

impl SlippageCurveConfig {
    /// The shape numbers the fit will try, evenly spaced over `[0, 1]` inclusive.
    ///
    /// **The grid lives here rather than inside the fit** so that a test asserting which shapes
    /// are reachable is asserting the same list the fit walks. One rung degenerates to the
    /// multiplying end alone, which is a configuration no caller should reach but which must not
    /// divide by zero.
    pub fn rise_shape_grid(&self) -> Vec<RiseShape> {
        let rungs = self.rise_shape_rungs.max(1);
        if rungs == 1 {
            return vec![RiseShape::MULTIPLYING];
        }
        (0..rungs)
            .map(|rung| {
                RiseShape::new(rung as f64 / (rungs - 1) as f64)
                    .expect("an evenly spaced rung of [0, 1] is a rise shape")
            })
            .collect()
    }
}

impl Default for SlippageCurveConfig {
    fn default() -> Self {
        Self {
            rise_shape_rungs: RISE_SHAPE_RUNGS,
            min_cells_for_a_curve: MIN_CELLS_FOR_A_CURVE,
            disagreement_knee: DISAGREEMENT_KNEE,
            draw_curves: true,
        }
    }
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_rise_shape_outside_zero_to_one_is_refused() {
        assert!(RiseShape::new(0.0).is_some());
        assert!(RiseShape::new(1.0).is_some());
        assert!(RiseShape::new(0.8).is_some());
        assert!(RiseShape::new(-0.01).is_none());
        assert!(RiseShape::new(1.01).is_none());
        assert!(RiseShape::new(f64::NAN).is_none());
        assert!(RiseShape::new(f64::INFINITY).is_none());
    }

    /// The two ends are the two arithmetics, and only exactly zero takes the logarithm branch.
    #[test]
    fn only_the_multiplying_end_takes_the_logarithm_branch() {
        assert!(RiseShape::MULTIPLYING.is_multiplying());
        assert!(!RiseShape::ADDING.is_multiplying());
        assert!(
            !RiseShape::new(0.05)
                .expect("inside the range")
                .is_multiplying()
        );
        assert_eq!(RiseShape::MULTIPLYING.get(), 0.0);
        assert_eq!(RiseShape::ADDING.get(), 1.0);
    }

    /// A cell's own precision is set by the reads that moved, not by the reads that crossed.
    #[test]
    fn a_cells_own_error_falls_as_the_square_root_of_its_slipped_reads() {
        let thin = FittedCell {
            repeats: 9,
            level: 0.0067,
            slipped_reads: 100.0,
        };
        let fat = FittedCell {
            repeats: 9,
            level: 0.0067,
            slipped_reads: 10_000.0,
        };
        assert!((thin.relative_standard_error() - 0.1).abs() < 1e-12);
        assert!((fat.relative_standard_error() - 0.01).abs() < 1e-12);
    }

    /// A cell the fit says nothing slipped in must not divide by zero.
    #[test]
    fn a_cell_with_no_slipped_reads_gets_an_error_of_one_rather_than_infinity() {
        let empty = FittedCell {
            repeats: 8,
            level: 0.0,
            slipped_reads: 0.0,
        };
        assert!(empty.relative_standard_error().is_finite());
        assert!((empty.relative_standard_error() - 1.0).abs() < 1e-12);
    }

    /// The three variants are one weight, so a consumer never matches on which end it is.
    #[test]
    fn every_level_source_reports_the_share_the_curve_carried() {
        assert_eq!(LevelSource::Cell.curve_weight(), 0.0);
        assert_eq!(LevelSource::Curve.curve_weight(), 1.0);
        assert_eq!(
            LevelSource::Blend { curve_weight: 0.93 }.curve_weight(),
            0.93
        );
    }

    /// A curve built at the adding end is a straight line in the level itself, and reading it
    /// back must return the line.
    #[test]
    fn the_adding_end_reads_back_as_a_straight_line_in_the_level() {
        // level = 0.005 * repeats - 0.035, which is 0.005 at 8 repeats and 0.115 at 30.
        let curve = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: -0.035,
            slope: 0.005,
            fitted_from: 8,
            fitted_to: 30,
            held_out_error: 0.044,
            cells: 23,
        };
        assert!((curve.level_at(8) - 0.005).abs() < 1e-12);
        assert!((curve.level_at(20) - 0.065).abs() < 1e-12);
        assert!((curve.level_at(30) - 0.115).abs() < 1e-12);
    }

    /// A curve built at the multiplying end multiplies by a fixed factor each repeat.
    #[test]
    fn the_multiplying_end_reads_back_as_a_fixed_factor_a_repeat() {
        // log(level) = -6.5 + 0.4 * repeats, so each repeat multiplies by exp(0.4) = 1.4918.
        let curve = SlippageCurve {
            rise_shape: RiseShape::MULTIPLYING,
            intercept: -6.5,
            slope: 0.4,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.124,
            cells: 5,
        };
        let at_eight = curve.level_at(8);
        let at_nine = curve.level_at(9);
        assert!(((-6.5_f64 + 0.4 * 8.0).exp() - at_eight).abs() < 1e-15);
        assert!((at_nine / at_eight - 0.4_f64.exp()).abs() < 1e-12);
    }

    /// Beyond the repeat counts the curve saw, the level is the nearer fitted end's — never the
    /// line continued. This is the property that stops an exponential reporting 205.
    #[test]
    fn outside_the_fitted_range_the_level_is_held_at_the_nearer_end() {
        let curve = SlippageCurve {
            rise_shape: RiseShape::MULTIPLYING,
            intercept: -6.5,
            slope: 0.4,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.124,
            cells: 5,
        };
        assert_eq!(curve.level_at(7), curve.level_at(8));
        assert_eq!(curve.level_at(1), curve.level_at(8));
        assert_eq!(curve.level_at(13), curve.level_at(12));
        assert_eq!(curve.level_at(150), curve.level_at(12));

        // Continued rather than held, this curve reaches exp(-6.5 + 0.4*30) = 244.7 at 30
        // repeats — the shape of the failure that made holding flat the rule — where held it
        // stays at its 12-repeat value of 0.1827.
        assert_eq!(curve.level_on_the_line(30.0), LEVEL_CEILING);
        assert!((curve.level_at(30) - (-6.5_f64 + 0.4 * 12.0).exp()).abs() < 1e-15);
        assert!((curve.level_at(30) - 0.1827).abs() < 1e-4);
    }

    /// The property every consumer of this curve relies on: a longer tract never comes back
    /// less slippery than a shorter one. It has to hold at both ends of the shape grid and in
    /// between, because the shape number is fitted and any rung can win.
    #[test]
    fn a_positive_slope_never_lets_a_longer_tract_come_back_less_slippery() {
        for shape in SlippageCurveConfig::default().rise_shape_grid() {
            let (low, high) = (0.004_f64, 0.12_f64);
            let (low_t, high_t) = if shape.is_multiplying() {
                (low.ln(), high.ln())
            } else {
                (low.powf(shape.get()), high.powf(shape.get()))
            };
            let slope = (high_t - low_t) / (30.0 - 8.0);
            let curve = SlippageCurve {
                rise_shape: shape,
                intercept: low_t - slope * 8.0,
                slope,
                fitted_from: 8,
                fitted_to: 30,
                held_out_error: 0.05,
                cells: 23,
            };
            let mut previous = 0.0;
            for repeats in 1..=60 {
                let level = curve.level_at(repeats);
                assert!(
                    level >= previous - 1e-15,
                    "shape {shape} fell from {previous} to {level} at {repeats} repeats"
                );
                assert!((LEVEL_FLOOR..=LEVEL_CEILING).contains(&level));
                previous = level;
            }
            assert!((curve.level_at(8) - low).abs() < 1e-12, "shape {shape}");
            assert!((curve.level_at(30) - high).abs() < 1e-12, "shape {shape}");
        }
    }

    #[test]
    #[should_panic(expected = "has no range to hold at its ends")]
    fn a_curve_whose_range_runs_backwards_says_so_rather_than_clamping() {
        let curve = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: 0.0,
            slope: 0.005,
            fitted_from: 30,
            fitted_to: 8,
            held_out_error: 0.044,
            cells: 23,
        };
        let _ = curve.level_at(12);
    }

    #[test]
    fn a_cell_knows_whether_it_sat_inside_the_curves_range() {
        let curve = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: -0.035,
            slope: 0.005,
            fitted_from: 8,
            fitted_to: 30,
            held_out_error: 0.044,
            cells: 23,
        };
        assert_eq!(curve.reach(7), CurveReach::BelowFitted);
        assert_eq!(curve.reach(8), CurveReach::Inside);
        assert_eq!(curve.reach(30), CurveReach::Inside);
        assert_eq!(curve.reach(31), CurveReach::AboveFitted);
    }

    /// A line fitted at the adding end goes negative below its range; the level may not.
    #[test]
    fn a_level_is_never_reported_outside_the_open_unit_interval() {
        let adding = SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: -0.035,
            slope: 0.005,
            fitted_from: 8,
            fitted_to: 30,
            held_out_error: 0.044,
            cells: 23,
        };
        // -0.035 + 0.005*5 = -0.010, a level a line is happy to produce and a caller cannot use.
        assert_eq!(adding.level_on_the_line(5.0), LEVEL_FLOOR);
        assert!(adding.level_on_the_line(1000.0) <= LEVEL_CEILING);

        let multiplying = SlippageCurve {
            rise_shape: RiseShape::MULTIPLYING,
            intercept: -6.5,
            slope: 0.4,
            fitted_from: 8,
            fitted_to: 12,
            held_out_error: 0.124,
            cells: 5,
        };
        assert_eq!(multiplying.level_on_the_line(100.0), LEVEL_CEILING);
        assert!(multiplying.level_on_the_line(-1000.0) >= LEVEL_FLOOR);
    }

    /// A shape number strictly between the two ends inverts through a root, and the round trip
    /// has to land back on the value the line was built from.
    #[test]
    fn a_shape_number_between_the_ends_round_trips_through_its_root() {
        let shape = RiseShape::new(0.8).expect("inside the range");
        // Build the line so that it passes exactly through 0.03 at 12 repeats and 0.06 at 25.
        let (low, high) = (0.03_f64.powf(0.8), 0.06_f64.powf(0.8));
        let slope = (high - low) / (25.0 - 12.0);
        let curve = SlippageCurve {
            rise_shape: shape,
            intercept: low - slope * 12.0,
            slope,
            fitted_from: 12,
            fitted_to: 25,
            held_out_error: 0.038,
            cells: 20,
        };
        assert!((curve.level_at(12) - 0.03).abs() < 1e-12);
        assert!((curve.level_at(25) - 0.06).abs() < 1e-12);
        // and it sits between the two ends: below a straight line, above an exponential.
        let midpoint = curve.level_at(18);
        assert!(midpoint > 0.03 && midpoint < 0.06);
    }

    /// Cells drawn exactly on a straight line in the level are recovered at the adding end.
    #[test]
    fn cells_lying_on_a_straight_line_are_recovered_at_the_adding_end() {
        let cells: Vec<FittedCell> = (8..=30)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.005 * repeats as f64 - 0.035,
                slipped_reads: 1_000.0,
            })
            .collect();
        let curve = fit_line(&cells, RiseShape::ADDING).expect("a rising line");
        assert!((curve.slope - 0.005).abs() < 1e-12, "slope {}", curve.slope);
        assert!((curve.intercept + 0.035).abs() < 1e-12);
        assert_eq!((curve.fitted_from, curve.fitted_to), (8, 30));
        assert_eq!(curve.cells, 23);
        for cell in &cells {
            assert!((curve.level_at(cell.repeats) - cell.level).abs() < 1e-12);
        }
    }

    /// Cells drawn exactly on an exponential are recovered at the multiplying end.
    #[test]
    fn cells_lying_on_an_exponential_are_recovered_at_the_multiplying_end() {
        let cells: Vec<FittedCell> = (8..=12)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.002 * 1.47_f64.powi(repeats as i32 - 8),
                slipped_reads: 500.0,
            })
            .collect();
        let curve = fit_line(&cells, RiseShape::MULTIPLYING).expect("a rising line");
        assert!((curve.slope - 1.47_f64.ln()).abs() < 1e-12);
        for cell in &cells {
            assert!((curve.level_at(cell.repeats) / cell.level - 1.0).abs() < 1e-12);
        }
    }

    /// A fit that says slippage falls with repeat count is refused, never emitted.
    #[test]
    fn a_falling_fit_is_refused_rather_than_handed_to_a_caller() {
        let falling: Vec<FittedCell> = (8..=12)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.05 - 0.003 * repeats as f64,
                slipped_reads: 1_000.0,
            })
            .collect();
        match fit_line(&falling, RiseShape::ADDING) {
            Err(NoCurve::LevelWouldFall { slope }) => assert!(slope < 0.0),
            other => panic!("a falling set gave {other:?}"),
        }
    }

    /// Every cell at one repeat count leaves no line to draw, and it is not a falling fit.
    #[test]
    fn cells_at_a_single_repeat_count_have_no_line_through_them() {
        let flat = vec![
            FittedCell {
                repeats: 9,
                level: 0.004,
                slipped_reads: 500.0,
            },
            FittedCell {
                repeats: 9,
                level: 0.006,
                slipped_reads: 500.0,
            },
        ];
        assert_eq!(
            fit_line(&flat, RiseShape::ADDING),
            Err(NoCurve::OneRepeatCountOnly { repeats: 9 })
        );
    }

    /// A cell whose level fitted to exactly zero has no logarithm and no evidence; it is dropped
    /// before the fit rather than crashing it.
    #[test]
    fn a_cell_at_a_level_of_zero_is_dropped_before_the_fit() {
        let mut cells: Vec<FittedCell> = (8..=12)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.002 * 1.47_f64.powi(repeats as i32 - 8),
                slipped_reads: 500.0,
            })
            .collect();
        let clean = fit_line(&cells, RiseShape::MULTIPLYING).expect("a rising line");
        cells.push(FittedCell {
            repeats: 13,
            level: 0.0,
            slipped_reads: 0.0,
        });
        let with_empty = fit_line(&cells, RiseShape::MULTIPLYING).expect("a rising line");
        assert_eq!(with_empty.cells, clean.cells);
        assert_eq!(with_empty.fitted_to, 12);
        assert!((with_empty.slope - clean.slope).abs() < 1e-15);
    }

    /// One usable cell is not a line, and the refusal names the arithmetic floor rather than the
    /// one that is about whether the answer means anything.
    #[test]
    fn one_usable_cell_is_refused_at_the_arithmetic_floor() {
        let single = vec![FittedCell {
            repeats: 9,
            level: 0.004,
            slipped_reads: 500.0,
        }];
        assert_eq!(
            fit_line(&single, RiseShape::ADDING),
            Err(NoCurve::TooFewCells { cells: 1, floor: 2 })
        );
    }

    /// The weight is the cell's slipped reads, so a cell with a thousand times the evidence
    /// pulls the line to itself and a near-empty cell barely moves it.
    #[test]
    fn a_cell_pulls_the_line_in_proportion_to_its_slipped_reads() {
        let on_the_line: Vec<FittedCell> = (8..=30)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.005 * repeats as f64 - 0.035,
                slipped_reads: 10_000.0,
            })
            .collect();
        let mut with_an_outlier = on_the_line.clone();
        with_an_outlier.push(FittedCell {
            repeats: 20,
            level: 0.5,
            slipped_reads: 1.0,
        });
        let clean = fit_line(&on_the_line, RiseShape::ADDING).expect("a rising line");
        let pulled = fit_line(&with_an_outlier, RiseShape::ADDING).expect("a rising line");
        let moved = (pulled.level_at(20) - clean.level_at(20)).abs() / clean.level_at(20);
        assert!(
            moved < 0.02,
            "one read moved the line at 20 repeats by {:.1}%",
            moved * 100.0
        );
    }

    /// Cells drawn on a straight line make the grid pick the adding end, and cells drawn on an
    /// exponential make it pick the multiplying end. Neither rung has more freedom than the
    /// other, so only held-out prediction can separate them.
    #[test]
    fn the_grid_finds_the_shape_the_cells_were_drawn_with() {
        let straight: Vec<FittedCell> = (8..=30)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.005 * repeats as f64 - 0.035,
                slipped_reads: 5_000.0,
            })
            .collect();
        let curves = choose_rise_shape(&[straight], &SlippageCurveConfig::default())
            .expect("twenty-three cells on a line");
        assert_eq!(curves.rise_shape, RiseShape::ADDING);
        assert!(curves.held_out_error < 1e-6, "{}", curves.held_out_error);

        let compounding: Vec<FittedCell> = (8..=30)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.002 * 1.15_f64.powi(repeats as i32 - 8),
                slipped_reads: 5_000.0,
            })
            .collect();
        let curves = choose_rise_shape(&[compounding], &SlippageCurveConfig::default())
            .expect("twenty-three cells on an exponential");
        assert_eq!(curves.rise_shape, RiseShape::MULTIPLYING);
        assert!(curves.held_out_error < 1e-6, "{}", curves.held_out_error);
    }

    /// The shape number is one number for the period; the lines are one per group. Two groups
    /// with genuinely different levels must keep their own lines and share the rung.
    #[test]
    fn two_groups_share_a_shape_number_and_keep_their_own_lines() {
        let group_of = |scale: f64| -> Vec<FittedCell> {
            (8..=30)
                .map(|repeats| FittedCell {
                    repeats,
                    level: scale * (0.005 * repeats as f64 - 0.035),
                    slipped_reads: 5_000.0,
                })
                .collect()
        };
        let curves = choose_rise_shape(
            &[group_of(1.0), group_of(2.0)],
            &SlippageCurveConfig::default(),
        )
        .expect("two groups of twenty-three cells");
        assert_eq!(curves.rise_shape, RiseShape::ADDING);
        let lines: Vec<SlippageCurve> = curves.by_group.iter().flatten().copied().collect();
        assert_eq!(lines.len(), 2);
        assert!((lines[1].slope / lines[0].slope - 2.0).abs() < 1e-9);
        assert_eq!(curves.cells, 46);
        // Both carry the period's own held-out error, so one curve can be weighed alone.
        assert_eq!(lines[0].held_out_error, curves.held_out_error);
        assert_eq!(lines[1].held_out_error, curves.held_out_error);
    }

    /// A group that put no cell in this period gets no line, and does not stop the others.
    #[test]
    fn a_group_with_no_cells_gets_no_line_and_does_not_block_the_period() {
        let cells: Vec<FittedCell> = (8..=30)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.005 * repeats as f64 - 0.035,
                slipped_reads: 5_000.0,
            })
            .collect();
        let curves = choose_rise_shape(
            &[Vec::new(), cells, Vec::new()],
            &SlippageCurveConfig::default(),
        )
        .expect("one group carries the period");
        assert!(curves.by_group[0].is_none());
        assert!(curves.by_group[1].is_some());
        assert!(curves.by_group[2].is_none());
    }

    /// Below the floor a period gets no curve at all, and the refusal names the floor that
    /// turned it away rather than the arithmetic one.
    #[test]
    fn a_period_below_the_cell_floor_is_refused_by_that_floor() {
        let three: Vec<FittedCell> = (8..=10)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.002 * repeats as f64,
                slipped_reads: 500.0,
            })
            .collect();
        assert_eq!(
            choose_rise_shape(&[three], &SlippageCurveConfig::default()),
            Err(NoCurve::TooFewCells {
                cells: 3,
                floor: MIN_CELLS_FOR_A_CURVE,
            })
        );
    }

    /// Cells whose level falls with repeat count give no curve at any rung.
    #[test]
    fn a_period_whose_cells_fall_gets_no_curve_at_any_rung() {
        let falling: Vec<FittedCell> = (8..=20)
            .map(|repeats| FittedCell {
                repeats,
                level: 0.15 - 0.005 * repeats as f64,
                slipped_reads: 5_000.0,
            })
            .collect();
        assert!(matches!(
            choose_rise_shape(&[falling], &SlippageCurveConfig::default()),
            Err(NoCurve::LevelWouldFall { .. })
        ));
    }

    /// Where two rungs score the same, the larger wins — the two ends fail differently, and the
    /// multiplying end is the one that returns numbers above one outside its range.
    #[test]
    fn an_exact_tie_between_rungs_goes_to_the_larger() {
        // Two cells lie exactly on every rung's curve, so every rung scores zero.
        let two_plus: Vec<FittedCell> = vec![
            FittedCell {
                repeats: 8,
                level: 0.004,
                slipped_reads: 1_000.0,
            },
            FittedCell {
                repeats: 12,
                level: 0.008,
                slipped_reads: 1_000.0,
            },
            FittedCell {
                repeats: 16,
                level: 0.012,
                slipped_reads: 1_000.0,
            },
            FittedCell {
                repeats: 20,
                level: 0.016,
                slipped_reads: 1_000.0,
            },
        ];
        let curves =
            choose_rise_shape(&[two_plus], &SlippageCurveConfig::default()).expect("four cells");
        // The straight line fits these exactly, so no rung can beat it and the tie rule keeps
        // the largest that ties.
        assert_eq!(curves.rise_shape, RiseShape::ADDING);
    }

    /// The curve HG002's homopolymers give, for the blend tests: a straight line whose held-out
    /// error is the 4.4% spec §7 quotes.
    fn hg002_homopolymer_curve() -> SlippageCurve {
        SlippageCurve {
            rise_shape: RiseShape::ADDING,
            intercept: -0.035,
            slope: 0.005,
            fitted_from: 8,
            fitted_to: 30,
            held_out_error: 0.044,
            cells: 23,
        }
    }

    /// **The two figures spec §7 quotes**, and the reason the blend is not a switch: the curve
    /// carries most of the weight at a near-empty cell and almost none at a full one.
    #[test]
    fn the_curve_carries_the_weight_at_a_thin_cell_and_stands_aside_at_a_full_one() {
        let curve = hg002_homopolymer_curve();
        let at = |slipped: f64| {
            let cell = FittedCell {
                repeats: 20,
                level: curve.level_at(20),
                slipped_reads: slipped,
            };
            blend_level(
                Some(cell),
                Some(&curve),
                20,
                &SlippageCurveConfig::default(),
            )
            .expect("a cell and a curve")
            .source
            .curve_weight()
        };
        let thin = at(40.0);
        let full = at(8_000.0);
        assert!(
            (thin - 0.93).abs() < 0.01,
            "at 40 slipped reads the curve carried {thin:.3}, and the spec quotes 0.93"
        );
        assert!(
            (full - 0.06).abs() < 0.01,
            "at 8,000 slipped reads the curve carried {full:.3}, and the spec quotes 0.06"
        );
        assert!(thin > full);
    }

    /// **The case that decides against using the curve everywhere.** HG002's 9-repeat
    /// homopolymer: the cell's own 3,520 slipped reads say 6.7 reads slipping per 1,000 and the
    /// curve says 10.4. The blend must stay near the cell, not near the curve.
    #[test]
    fn a_well_measured_cell_keeps_its_answer_where_the_curve_is_wrong_about_it() {
        let curve = hg002_homopolymer_curve();
        let cell = FittedCell {
            repeats: 9,
            level: 0.00673,
            slipped_reads: 3_520.0,
        };
        assert!(
            (curve.level_at(9) - 0.0104).abs() < 5e-4,
            "the curve says {} at 9 repeats",
            curve.level_at(9)
        );

        let with_knee =
            blend_level(Some(cell), Some(&curve), 9, &SlippageCurveConfig::default()).unwrap();
        let no_knee = blend_level(
            Some(cell),
            Some(&curve),
            9,
            &SlippageCurveConfig {
                disagreement_knee: f64::INFINITY,
                ..SlippageCurveConfig::default()
            },
        )
        .unwrap();

        // Always-curve would report 10.4 per 1,000; both blends stay within a tenth of that of
        // the cell's own 6.73.
        assert!((no_knee.level - 0.00712).abs() < 2e-4, "{}", no_knee.level);
        assert!(
            (with_knee.level - 0.00676).abs() < 2e-4,
            "{}",
            with_knee.level
        );
        assert!(
            with_knee.level < no_knee.level,
            "the knee should pull further toward the cell"
        );
    }

    /// The knee only fires on a gap no sampling noise produces; a cell sitting close to its
    /// curve is blended as if the knee were not there.
    #[test]
    fn the_knee_leaves_a_cell_that_agrees_with_its_curve_alone() {
        let curve = hg002_homopolymer_curve();
        let cell = FittedCell {
            repeats: 20,
            // 3% off the curve, which at these errors is well inside the knee.
            level: curve.level_at(20) * 1.03,
            slipped_reads: 3_671.0,
        };
        let with_knee = blend_level(
            Some(cell),
            Some(&curve),
            20,
            &SlippageCurveConfig::default(),
        )
        .unwrap();
        let no_knee = blend_level(
            Some(cell),
            Some(&curve),
            20,
            &SlippageCurveConfig {
                disagreement_knee: f64::INFINITY,
                ..SlippageCurveConfig::default()
            },
        )
        .unwrap();
        assert_eq!(with_knee, no_knee);
    }

    /// The three degenerate cases, none of which needs a branch at the call site.
    #[test]
    fn the_three_degenerate_cases_fall_out_of_the_same_formula() {
        let curve = hg002_homopolymer_curve();
        let cell = FittedCell {
            repeats: 20,
            level: 0.0675,
            slipped_reads: 3_671.0,
        };

        // A cell with no curve keeps its own level, exactly as today.
        let alone = blend_level(Some(cell), None, 20, &SlippageCurveConfig::default()).unwrap();
        assert_eq!(alone.source, LevelSource::Cell);
        assert_eq!(alone.level, cell.level);
        assert!(alone.reach.is_none());

        // A cell with no fit of its own takes the curve whole — the case this change exists for.
        let borrowed =
            blend_level(None, Some(&curve), 14, &SlippageCurveConfig::default()).unwrap();
        assert_eq!(borrowed.source, LevelSource::Curve);
        assert_eq!(borrowed.level, curve.level_at(14));
        assert_eq!(borrowed.reach, Some(CurveReach::Inside));

        // Neither is nothing, and it says so rather than emitting a zero.
        assert!(blend_level(None, None, 14, &SlippageCurveConfig::default()).is_none());
    }

    /// A cell the fit put at exactly zero has no logarithm; it takes the curve rather than
    /// dragging the blend to zero.
    #[test]
    fn a_cell_fitted_at_zero_takes_the_curve_whole() {
        let curve = hg002_homopolymer_curve();
        let empty = FittedCell {
            repeats: 12,
            level: 0.0,
            slipped_reads: 0.0,
        };
        let blended = blend_level(
            Some(empty),
            Some(&curve),
            12,
            &SlippageCurveConfig::default(),
        )
        .unwrap();
        assert_eq!(blended.source, LevelSource::Curve);
        assert_eq!(blended.level, curve.level_at(12));
    }

    /// A cell above the curve's fitted range takes a held level and is told so.
    #[test]
    fn a_cell_beyond_the_curves_range_is_marked_as_held() {
        let curve = hg002_homopolymer_curve();
        let blended = blend_level(None, Some(&curve), 45, &SlippageCurveConfig::default()).unwrap();
        assert_eq!(blended.reach, Some(CurveReach::AboveFitted));
        assert_eq!(blended.level, curve.level_at(30));
    }

    /// Whatever the inputs, an emitted level is a probability.
    #[test]
    fn an_emitted_level_is_always_a_probability() {
        let curve = hg002_homopolymer_curve();
        for slipped in [0.5, 1.0, 40.0, 3_500.0, 100_000.0] {
            for level in [1e-9, 1e-4, 0.02, 0.5, 0.98] {
                for repeats in [1_u64, 9, 20, 30, 120] {
                    let cell = FittedCell {
                        repeats,
                        level,
                        slipped_reads: slipped,
                    };
                    let blended = blend_level(
                        Some(cell),
                        Some(&curve),
                        repeats,
                        &SlippageCurveConfig::default(),
                    )
                    .expect("a cell and a curve");
                    assert!(
                        (LEVEL_FLOOR..=LEVEL_CEILING).contains(&blended.level),
                        "level {level} at {repeats} repeats with {slipped} slipped reads gave {}",
                        blended.level
                    );
                    assert!((0.0..=1.0).contains(&blended.source.curve_weight()));
                }
            }
        }
    }

    #[test]
    fn the_default_config_is_the_measured_constants_with_curves_on() {
        let config = SlippageCurveConfig::default();
        assert_eq!(config.rise_shape_rungs, RISE_SHAPE_RUNGS);
        assert_eq!(config.min_cells_for_a_curve, MIN_CELLS_FOR_A_CURVE);
        assert_eq!(config.disagreement_knee, DISAGREEMENT_KNEE);
        assert!(config.draw_curves);
    }

    /// The grid must be able to name the answers both cohorts gave, or the fit cannot return
    /// them: 21 rungs over `[0, 1]` are steps of 0.05, and 0.80 is one of them.
    #[test]
    fn the_shape_grid_resolves_the_answers_both_cohorts_gave() {
        let grid = SlippageCurveConfig::default().rise_shape_grid();
        assert_eq!(grid.len(), RISE_SHAPE_RUNGS);
        for wanted in [0.0, 0.8, 1.0] {
            assert!(
                grid.iter().any(|shape| (shape.get() - wanted).abs() < 1e-9),
                "the grid cannot name {wanted}"
            );
        }
        assert_eq!(grid.first().copied(), Some(RiseShape::MULTIPLYING));
        assert_eq!(grid.last().copied(), Some(RiseShape::ADDING));
    }

    /// A grid of one rung must not divide by zero when spacing the rungs.
    #[test]
    fn a_single_rung_grid_is_the_multiplying_end_alone() {
        let config = SlippageCurveConfig {
            rise_shape_rungs: 1,
            ..SlippageCurveConfig::default()
        };
        assert_eq!(config.rise_shape_grid(), vec![RiseShape::MULTIPLYING]);
        let none = SlippageCurveConfig {
            rise_shape_rungs: 0,
            ..SlippageCurveConfig::default()
        };
        assert_eq!(none.rise_shape_grid(), vec![RiseShape::MULTIPLYING]);
    }
}
