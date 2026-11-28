//! Where an emitted slippage number came from: a stratum's own fit, or a point on its motif
//! length's penalised curve.
//!
//! **The three slippage numbers.** When a read crosses a repeat tract, the copying steps before
//! sequencing sometimes add or drop whole repeat units, so the read reports a tract of the wrong
//! length. The *level* is how often a read slips at all; of the reads that slipped, the
//! **direction split** is the share that came back *shorter* rather than longer, and the
//! **fall-off** is how much rarer a two-unit slip is than a one-unit slip.
//!
//! [`LevelSource`] says where a level came from and [`ShareSource`] where one of the two shares
//! did; a point on a penalised curve carries its [`PenalisedPoint`], and [`CurveReach`] records
//! whether the stratum sat inside the range the curve was fitted over.
//!
//! **Only these sources are left.** The repeat-tract fit before the penalised one (the climb,
//! deleted on 2026-10-16) also read numbers off a parametric curve drawn across repeat count, or
//! blended a stratum's own number with one. A parameters file that records either is refused
//! when it is read (`ParametersFile::from_toml`).

// ---------------------------------------------------------------------
// Where an emitted number came from
// ---------------------------------------------------------------------

/// Where a cell's emitted slippage level came from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LevelSource {
    /// The cell's own fit, whole.
    Cell,
    /// **A point on its motif length's penalised curve, fitted with this cell's own reads.**
    ///
    /// The penalised fit fits every repeat count of a motif length together, pulling each
    /// count's level toward its neighbours'. The fit reports only cells with reads, so the value
    /// is never read across a gap; it is a variant of its own rather than [`LevelSource::Cell`]
    /// because its neighbours' reads pulled it too (`str_fit_to_calling.md` §2.3).
    PenalisedCurve(PenalisedPoint),
}

/// Where one of a stratum's two emitted shares came from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ShareSource {
    /// The stratum's own fit, whole.
    Stratum,
    /// **A point on its motif length's penalised curve, fitted with this stratum's own reads** —
    /// as for [`LevelSource::PenalisedCurve`].
    PenalisedCurve(PenalisedPoint),
}

/// **What the parameters file records about one point on a penalised curve.**
///
/// The penalised fit draws one curve per motif length through every repeat count with reads; how
/// stiff the curve is was chosen from the data, and at some repeat counts the fit let the curve
/// bend more sharply because the reads insisted on a kink. Both facts travel with each number so
/// that a reader of the file can tell a value held close to its neighbours from one the reads set
/// on their own.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PenalisedPoint {
    /// The curve's stiffness strength: larger holds each repeat count's value closer to its
    /// neighbours'.
    pub strength: f64,
    /// Whether a bend was eased at this repeat count: the reads insisted on a kink here.
    pub eased_here: bool,
}

/// Whether a stratum sat inside the repeat counts its curve was fitted over.
///
/// The penalised fit reports only strata with reads, so it always writes [`CurveReach::Inside`];
/// the other two are the file's vocabulary, which it still reads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CurveReach {
    /// Inside the fitted range.
    Inside,
    /// Below the lowest stratum that fed the curve.
    BelowFitted,
    /// Above the highest.
    AboveFitted,
}
