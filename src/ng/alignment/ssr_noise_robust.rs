//! Algorithm 4n — the **noise-robust** repeat delimiter.
//!
//! Same skeleton as algorithm 4 ([`super::ssr_best_path_unit_slip`]): a five-state Viterbi
//! (match / out-of-frame insertion / out-of-frame deletion / whole-unit slip insertion /
//! whole-unit slip deletion) whose slip costs come from the shared [`StutterModel`]. What
//! differs is **where the tract boundary is allowed to move**, and the reason is what the
//! synthetic bake-off measures: a couple of miscalled bases near a flank/tract junction are
//! enough to make a *shifted* line-up score best, and the delimiter then reports a length that
//! is off by a base or a whole unit.
//!
//! The arithmetic behind that is worth stating plainly, because it is the whole design. At
//! Q40 a mismatch scores `ln(1e-4/3) ≈ -10.31`, and algorithm 4's flank gap-open scores
//! `ln(2.9e-5) ≈ -10.45`. **One bad base costs the same as opening a gap.** So an insertion
//! and a deletion — which together slide the read one base against the reference — cost about
//! as much as two mismatches, and any noisy read where a one-base slide repairs three of them
//! prefers the slide. The slide moves the junction, and the measurement is wrong by a base.
//!
//! [`NoiseRobustConfig`] holds the four knobs that answer this, all **constructor state**, so
//! which regime an experiment ran at is visible at the call site rather than hidden in a
//! scratch buffer:
//!
//! 1. **`flank_gaps`** — whether an ordinary per-base gap may open outside the tract. Off by
//!    default. A read's flank has no reason to carry an indel: the flank *is* the anchor, and
//!    its job is to hold the boundary still. Banning gaps there prices a bad flank base as
//!    what it is — one mismatch — and means **no single base, and no pair of bases, can move
//!    the junction**. That is the "anchoring supported by more than one base" idea in its
//!    strongest form: the entire flank anchors, jointly.
//!
//! 2. **`junction_guard_bases`** — how many tract columns at each end may not *open* an
//!    out-of-frame gap. Six by default, one motif unit at the widest period, and always
//!    clamped so the tract keeps an interior. The out-of-frame route exists for genuinely
//!    impure tracts; near the junction it does something else, trading one mismatch for a
//!    one-base boundary shift. Whole-unit slips remain free to open anywhere in the tract, so
//!    real length variation is untouched — only the sub-unit nudge at the boundary goes.
//!
//! 3. **`min_flank_support`** — how many read bases must be aligned against a flank before
//!    that flank counts as anchoring. Algorithm 4 accepts a junction crossed by a *single*
//!    base, so a read that ran off deep inside a long tract can pick up one stray base on the
//!    far flank and be reported as a confident short measurement — the exact failure
//!    [`RepeatSpan`]'s partial variants exist to prevent. This is the direct form of the same
//!    idea as knob 1: there, no one base may *move* the boundary; here, no one base may
//!    *declare* it.
//!
//! 4. **`slip_margin`** — an extra log-cost charged once per whole-unit slip run. Under heavy
//!    miscalling a spurious slip can win, because the mismatching bases it swallows get
//!    re-scored against the motif instead of against the reference. A fixed margin leaves
//!    genuine stutter alone (it recovers a whole unit of matched bases, far more than the
//!    margin) while making a noise-driven slip lose. It ships at zero — see
//!    [`NoiseRobustConfig::default`].
//!
//! **Two routes are deliberately exempt from the flank ban**, and getting this wrong would
//! destroy the partial cases rather than the noisy ones. A read that stops inside the tract
//! still has to reach the frame's end, and one that starts inside the tract still has to be
//! reached from the frame's start; both cross a flank by deletion. Those crossings are exactly
//! what [`RepeatSpan::FromLeft`] and [`RepeatSpan::FromRight`] report, so they keep the
//! ordinary open cost. They are recognised structurally — a deletion on read row 0 (nothing
//! consumed yet) or on the last read row (nothing left to consume) — never by a threshold.
//!
//! Everything else is algorithm 4's, unchanged: the recurrence, the emission model, the
//! tie-break order, the traceback and the [`RepeatSpan`] classification. The two are therefore
//! directly comparable, and error-free reads measure identically under both.

use super::emission::Emission;
use super::ssr_best_path_flat_gap::{TractReadout, TransitionCosts};
use super::stutter::StutterModel;
use super::{BestPathAligner, ReadBases, RepeatContext, RepeatSpan};

/// The five states a cell can be entered in — algorithm 4's, unchanged.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum State {
    /// A read base and a reference base consumed together.
    Match = 0,
    /// One read base against no reference base — an out-of-frame insertion.
    Insertion = 1,
    /// One reference base against no read base — an out-of-frame deletion.
    Deletion = 2,
    /// A whole unit of read bases the reference lacks — an in-frame expansion.
    SlipInsertion = 3,
    /// A whole unit of reference bases the read lacks — an in-frame contraction.
    SlipDeletion = 4,
}

const STATES: usize = 5;

impl State {
    #[inline]
    const fn index(self) -> usize {
        self as usize
    }
}

/// An unreachable score. A real value: some cells cannot be entered in some states, and this
/// aligner also uses it to *close* a route (a banned gap-open is an unreachable transition,
/// not an expensive one).
const UNREACHABLE: f64 = f64::NEG_INFINITY;

/// The four knobs that separate this delimiter from algorithm 4 — see the module docs for
/// what each one is for.
///
/// It is **constructor state, not scratch**: each value changes results, so a bake-off has to
/// be able to see and report the setting it compared at. [`Default`] is the shipped
/// configuration, and [`SsrNoiseRobustAligner::new`] uses it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NoiseRobustConfig {
    /// Whether an ordinary per-base gap may open **outside** the tract. `false` ships.
    pub flank_gaps: bool,
    /// How many tract columns at each end are barred from *opening* an out-of-frame gap.
    /// `6` ships — one motif unit at the widest period [`Motif`] can hold.
    ///
    /// A width rather than a flag because the right number is an empirical question, and a
    /// reader who sees a `bool` cannot tell whether widening it was tried. It was, and the
    /// synthetic bake-off keeps rewarding a wider guard right up to the point where the tract
    /// has no interior left — which is the tell that past a unit or so it is no longer guarding
    /// a junction but simply deleting the out-of-frame route. Six bases is where the
    /// justification stops (a sub-unit gap within one unit of the junction is a boundary
    /// convention, not an interior event), so that is where the default stops.
    ///
    /// **Always clamped so the tract keeps an interior**: see [`effective_guard`], which leaves
    /// at least one motif unit of openable columns in the middle. Without it a short tract —
    /// three or four units, the common case in a real genome — would have every out-of-frame
    /// route closed, and an interrupted repeat there would be forced onto a tidy unit count.
    ///
    /// [`Motif`]: crate::ng::types::Motif
    pub junction_guard_bases: usize,
    /// Extra log-cost per whole-unit slip run, in nats. `0.0` ships — see
    /// [`Self::default`] for why.
    pub slip_margin: f64,
    /// How many read bases must be **aligned against a flank** before that flank counts as
    /// anchoring the tract. `3` ships.
    ///
    /// Algorithm 4 treats a flank as anchoring whenever the traceback crossed its junction at
    /// all, so a *single* read base landing on the first flank column is enough to turn a read
    /// that ran off inside the tract into a confident measurement — a long allele reported as a
    /// short one, which is the worst answer the span type exists to prevent. Requiring several
    /// bases is the counterpart to banning flank gaps: there, no one base may *move* the
    /// boundary; here, no one base may *declare* it.
    ///
    /// Counted in aligned columns, not in agreeing bases, so the rule is unaffected by
    /// miscalls — a noisy flank still anchors, a one-base flank still does not. The value is
    /// clamped to the flank that actually exists, so a repeat near a contig edge is judged
    /// against the flank it has rather than being made unanchorable.
    ///
    /// **Worth being straight about what the measurement says.** On the synthetic bake-off
    /// this knob is decisive at a narrow junction guard (it is what takes `partial` from
    /// 0.905 to 1.000 at `junction_guard_bases: 1`) but makes **no difference at all** at the
    /// shipped guard of six, which already closes the same door from the other side. It ships
    /// on regardless, because the two are not the same guarantee: the guard makes a one-base
    /// anchor *unlikely to win a scoring contest*, while this makes it *inadmissible*. The
    /// failure it blocks — a long allele reported as a confident short measurement — is the
    /// one wrong answer in this module that a downstream genotype cannot recover from, and a
    /// belt that costs one comparison per read is worth wearing next to the braces.
    ///
    /// Three, not five, so it stays below production's own five-base flank admission
    /// threshold: this rule is meant to reject the degenerate one- and two-base anchors, not
    /// to relitigate where a read is admitted.
    pub min_flank_support: usize,
}

impl Default for NoiseRobustConfig {
    /// The shipped configuration, chosen on the synthetic bake-off.
    ///
    /// `slip_margin` is zero because the measurement did not support it: the slip route is
    /// already priced by the stutter model, and taxing it further cost more on genuine
    /// expansions than it recovered on noise. The knob stays because it is the natural third
    /// lever and the next data set may disagree — a zero here is a *finding*, not an absent
    /// feature.
    fn default() -> Self {
        Self {
            flank_gaps: false,
            junction_guard_bases: 6,
            slip_margin: 0.0,
            min_flank_support: 3,
        }
    }
}

/// How wide the junction guard may actually be at this locus: the configured width, cut back
/// so the tract keeps **at least one motif unit** of columns that can still open an
/// out-of-frame gap.
///
/// The clamp is the difference between a guard and a ban. A tract of three or four units — the
/// common case in a real genome, and shorter than anything the synthetic set exercises — is
/// narrower than two guards, so an unclamped six-base guard would close its whole interior and
/// force every interrupted repeat there onto a tidy unit count. That is precisely the failure
/// algorithm 4's out-of-frame route exists to avoid, so it must not be reintroduced as a side
/// effect of a boundary fix.
#[inline]
fn effective_guard(configured: usize, tract_len: usize, period: usize) -> usize {
    configured.min(tract_len.saturating_sub(period) / 2)
}

/// The whole-unit slip transition costs, in log space, derived from the shared
/// [`StutterModel`] — algorithm 4's derivation plus [`NoiseRobustConfig::slip_margin`].
#[derive(Debug, Clone, Copy)]
struct SlipCosts {
    /// `ln(in_up · in_geom) − ln(equal) − margin`.
    open_expansion: f64,
    /// `ln(in_down · in_geom) − ln(equal) − margin`.
    open_contraction: f64,
    /// `ln(1 − in_geom)`, per unit after the first.
    extend: f64,
}

impl SlipCosts {
    fn from_model(model: &StutterModel, margin: f64) -> Self {
        let ln_equal = model.equal().ln();
        Self {
            open_expansion: (model.in_up() * model.in_geom()).ln() - ln_equal - margin,
            open_contraction: (model.in_down() * model.in_geom()).ln() - ln_equal - margin,
            extend: (1.0 - model.in_geom()).ln(),
        }
    }
}

/// Per-worker scratch: a ring of `period + 1` score rows plus the full backpointer matrix,
/// exactly as algorithm 4 needs (a whole-unit insertion reaches `(i, j)` from `(i − period, j)`).
///
/// Grow-and-keep; buffers only, deciding nothing that changes a result.
#[derive(Debug, Default)]
pub struct NoiseRobustScratch {
    rows: Vec<Vec<[f64; STATES]>>,
    backpointers: Vec<[State; STATES]>,
}

impl NoiseRobustScratch {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    fn resize(&mut self, read_len: usize, reference_len: usize, ring_len: usize) {
        let width = reference_len + 1;
        if self.rows.len() < ring_len {
            self.rows.resize_with(ring_len, Vec::new);
        }
        for row in &mut self.rows {
            if row.len() < width {
                row.resize(width, [UNREACHABLE; STATES]);
            }
        }
        let cells = (read_len + 1) * width;
        if self.backpointers.len() < cells {
            self.backpointers.resize(cells, [State::Match; STATES]);
        }
    }
}

/// Pick the best candidate, keeping the **first on ties** — the spec's determinism rule,
/// encoded by the caller's argument order (match, then deletion, then insertion, then slips).
#[inline]
fn best_of(candidates: &[(f64, State)]) -> (f64, State) {
    debug_assert!(!candidates.is_empty(), "best_of requires a candidate");
    let mut best = candidates[0];
    for &candidate in &candidates[1..] {
        if candidate.0 > best.0 {
            best = candidate;
        }
    }
    best
}

/// Algorithm 4n: the noise-robust repeat delimiter. Generic over its [`Emission`], never
/// behind `dyn` (arch §4).
#[derive(Debug, Clone, Copy)]
pub struct SsrNoiseRobustAligner<E: Emission> {
    emission: E,
    costs: TransitionCosts,
    config: NoiseRobustConfig,
}

impl<E: Emission> SsrNoiseRobustAligner<E> {
    /// Build the delimiter at its shipped configuration.
    #[must_use]
    pub fn new(emission: E) -> Self {
        Self::with_config(emission, NoiseRobustConfig::default())
    }

    /// Build the delimiter at an explicit configuration — the form a bake-off uses, so the
    /// setting under comparison is named at the call site.
    #[must_use]
    pub fn with_config(emission: E, config: NoiseRobustConfig) -> Self {
        Self {
            emission,
            costs: TransitionCosts::new(),
            config,
        }
    }

    /// The configuration this delimiter is running at.
    #[must_use]
    pub fn config(&self) -> NoiseRobustConfig {
        self.config
    }

    /// Measure the read's repeat. `None` when there is no reference frame — a defined answer,
    /// not an error (arch §3).
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn delimit(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: &RepeatContext<'_>,
        scratch: &mut NoiseRobustScratch,
    ) -> Option<TractReadout> {
        debug_assert!(
            context
                .geometry
                .fits_reference(crate::ng::types::Bp(reference.len() as u64)),
            "the repeat geometry does not fit the reference stretch"
        );

        let read_len = read.len();
        let bases = read.bases();
        let reference_len = reference.len();
        if reference_len == 0 {
            return None;
        }

        let left_flank_len = usize::try_from(context.geometry.left_flank_len.get())
            .unwrap_or(usize::MAX)
            .min(reference_len);
        let right_flank_len = usize::try_from(context.geometry.right_flank_len.get())
            .unwrap_or(usize::MAX)
            .min(reference_len - left_flank_len);
        let tract_last = reference_len - right_flank_len;

        let motif = context.geometry.motif.as_bytes();
        let period = motif.len();
        debug_assert_eq!(period, context.geometry.motif.period());

        let slip = SlipCosts::from_model(context.stutter, self.config.slip_margin);
        let ring_len = period + 1;
        scratch.resize(read_len, reference_len, ring_len);
        let stride = reference_len + 1;
        let insertion_emission = self.emission.insert_ln();

        // A reference column is inside the tract when a gap touching it is a tract gap.
        let column_in_tract = |column: usize| left_flank_len < column && column <= tract_last;

        // The per-base gap-open cost by column, under the two boundary knobs: the tract-grade
        // open inside the tract, `UNREACHABLE` in the flank when `flank_gaps` is off, and
        // `UNREACHABLE` at the two junction columns when `junction_guard` is on.
        let guard = effective_guard(
            self.config.junction_guard_bases,
            tract_last - left_flank_len,
            period,
        );
        let gap_open = |column: usize| {
            if column_in_tract(column) {
                let near_junction = column <= left_flank_len + guard || column + guard > tract_last;
                if guard > 0 && near_junction {
                    UNREACHABLE
                } else {
                    self.costs.ln_gap_open_tract()
                }
            } else if self.config.flank_gaps {
                self.costs.ln_gap_open()
            } else {
                UNREACHABLE
            }
        };

        // The **terminal** deletion route, exempt from the flank ban: a read that stops inside
        // the tract still has to reach the frame's end, and one that starts inside it still has
        // to be reached from the frame's start. Those are the `FromLeft`/`FromRight` cases, and
        // closing them would turn every partial read into a wrong answer.
        let terminal_gap_open = self.costs.ln_gap_open();

        // The `k`-th base of a whole unit inserted at tract column `column`: the motif,
        // phase-aligned so a unit inserted at a unit boundary begins at `motif[0]`.
        let motif_base = |column: usize, k: usize| motif[((column - left_flank_len) + k) % period];

        let NoiseRobustScratch { rows, backpointers } = scratch;

        // Row 0 — no read base consumed. Every deletion here is a *leading* one (the read starts
        // later in the frame), so it takes the terminal cost regardless of the flank ban.
        {
            let row0 = &mut rows[0];
            row0[0] = [0.0, UNREACHABLE, UNREACHABLE, UNREACHABLE, UNREACHABLE];
            backpointers[0] = [State::Match; STATES];
            for column in 1..=reference_len {
                let (d, d_pred) = best_of(&[
                    (
                        terminal_gap_open + row0[column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend() + row0[column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                ]);
                let (sd, sd_pred) =
                    if column >= period && (column - period + 1..=column).all(column_in_tract) {
                        best_of(&[
                            (
                                slip.open_contraction + row0[column - period][State::Match.index()],
                                State::Match,
                            ),
                            (
                                slip.extend + row0[column - period][State::SlipDeletion.index()],
                                State::SlipDeletion,
                            ),
                        ])
                    } else {
                        (UNREACHABLE, State::Match)
                    };
                row0[column] = [UNREACHABLE, UNREACHABLE, d, UNREACHABLE, sd];
                backpointers[column] = [State::Match, State::Match, d_pred, State::Match, sd_pred];
            }
        }

        for row_index in 1..=read_len {
            let read_base = bases[row_index - 1];
            let scores = self.emission.scores_for(read.quality_at(row_index - 1));
            let back_row = row_index * stride;

            let prev_slot = (row_index - 1) % ring_len;
            let slip_slot = row_index.checked_sub(period).map(|r| r % ring_len);
            let cur_slot = row_index % ring_len;
            // A deletion on the last read row is a *trailing* one — the read has run out and the
            // rest of the frame must still be crossed. That is the `FromLeft` case, not the
            // boundary-sliding gap the flank ban is aimed at.
            let trailing = row_index == read_len;

            // Column 0 — a read base before any reference base: the leading-insertion route,
            // which is how a read that overhangs the frame's start is placed. Kept at the
            // ordinary open cost for the same reason the terminal deletions are.
            let (ins0, ins0_pred) = best_of(&[
                (
                    self.costs.ln_gap_open() + rows[prev_slot][0][State::Match.index()],
                    State::Match,
                ),
                (
                    self.costs.ln_gap_extend() + rows[prev_slot][0][State::Insertion.index()],
                    State::Insertion,
                ),
            ]);
            rows[cur_slot][0] = [
                UNREACHABLE,
                insertion_emission + ins0,
                UNREACHABLE,
                UNREACHABLE,
                UNREACHABLE,
            ];
            backpointers[back_row] = [
                State::Match,
                ins0_pred,
                State::Match,
                State::Match,
                State::Match,
            ];

            for column in 1..=reference_len {
                let emit = scores.pick(read_base, reference[column - 1]);

                // Match: from the diagonal, any state (a match may follow a completed slip).
                let (m, m_pred) = best_of(&[
                    (
                        self.costs.ln_match_to_match()
                            + rows[prev_slot][column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_close()
                            + rows[prev_slot][column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                    (
                        self.costs.ln_gap_close()
                            + rows[prev_slot][column - 1][State::Insertion.index()],
                        State::Insertion,
                    ),
                    (
                        rows[prev_slot][column - 1][State::SlipInsertion.index()],
                        State::SlipInsertion,
                    ),
                    (
                        rows[prev_slot][column - 1][State::SlipDeletion.index()],
                        State::SlipDeletion,
                    ),
                ]);

                let open = gap_open(column);

                // Out-of-frame single-base insertion (a read base with no reference base).
                let (ins, ins_pred) = best_of(&[
                    (
                        open + rows[prev_slot][column][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend()
                            + rows[prev_slot][column][State::Insertion.index()],
                        State::Insertion,
                    ),
                ]);

                // Out-of-frame single-base deletion (a reference base with no read base). On the
                // last read row the terminal route applies, so a banned flank open reopens.
                let del_open = if trailing {
                    terminal_gap_open.max(open)
                } else {
                    open
                };
                let (del, del_pred) = best_of(&[
                    (
                        del_open + rows[cur_slot][column - 1][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend()
                            + rows[cur_slot][column - 1][State::Deletion.index()],
                        State::Deletion,
                    ),
                ]);

                // Whole-unit insertion: `period` read bases forming a unit, scored against the
                // motif **each at its own quality**, at a tract column. Reachable from row
                // `i − period`, same column.
                let (sins, sins_pred) = if let Some(slip_slot) = slip_slot {
                    if column_in_tract(column) {
                        let unit_emit: f64 = (0..period)
                            .map(|k| {
                                let idx = row_index - period + k;
                                self.emission
                                    .scores_for(read.quality_at(idx))
                                    .pick(bases[idx], motif_base(column, k))
                            })
                            .sum();
                        let (score, pred) = best_of(&[
                            (
                                slip.open_expansion + rows[slip_slot][column][State::Match.index()],
                                State::Match,
                            ),
                            (
                                slip.extend + rows[slip_slot][column][State::SlipInsertion.index()],
                                State::SlipInsertion,
                            ),
                        ]);
                        (unit_emit + score, pred)
                    } else {
                        (UNREACHABLE, State::Match)
                    }
                } else {
                    (UNREACHABLE, State::Match)
                };

                // Whole-unit deletion: `period` reference bases the read lacks, when the deleted
                // unit lies wholly in the tract. No emissions — no read base is consumed.
                let (sdel, sdel_pred) =
                    if column >= period && (column - period + 1..=column).all(column_in_tract) {
                        best_of(&[
                            (
                                slip.open_contraction
                                    + rows[cur_slot][column - period][State::Match.index()],
                                State::Match,
                            ),
                            (
                                slip.extend
                                    + rows[cur_slot][column - period][State::SlipDeletion.index()],
                                State::SlipDeletion,
                            ),
                        ])
                    } else {
                        (UNREACHABLE, State::Match)
                    };

                rows[cur_slot][column] = [emit + m, insertion_emission + ins, del, sins, sdel];
                backpointers[back_row + column] =
                    [m_pred, ins_pred, del_pred, sins_pred, sdel_pred];
            }
        }

        // Final cell (m, n): best terminal state. Tie-break M > D > I > SlipIns > SlipDel.
        let last = &rows[read_len % ring_len][reference_len];
        let (_, final_state) = best_of(&[
            (last[State::Match.index()], State::Match),
            (last[State::Deletion.index()], State::Deletion),
            (last[State::Insertion.index()], State::Insertion),
            (last[State::SlipInsertion.index()], State::SlipInsertion),
            (last[State::SlipDeletion.index()], State::SlipDeletion),
        ]);

        let geometry = MatrixGeometry {
            stride,
            read_len,
            reference_len,
            left_flank_len,
            right_flank_len,
            period,
            min_flank_support: self.config.min_flank_support,
        };
        Some(trace_back(final_state, backpointers, geometry))
    }
}

/// The matrix shape the traceback needs — bundled so the walk takes one `Copy` argument rather
/// than six positional `usize`s, several of them interchangeable.
#[derive(Debug, Clone, Copy)]
struct MatrixGeometry {
    stride: usize,
    read_len: usize,
    reference_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
    period: usize,
    /// [`NoiseRobustConfig::min_flank_support`], carried into the walk that applies it.
    min_flank_support: usize,
}

/// Walk the traceback, read the tract off the two flank junctions, and **check that each flank
/// actually held it** — algorithm 4's walk plus the minimum-support rule.
///
/// The support counts are read bases *aligned against* a flank column, not bases that agree
/// with it: a miscalled flank base still holds the anchor, which is the whole point, while a
/// flank represented by one lone read base does not.
fn trace_back(
    final_state: State,
    backpointers: &[[State; STATES]],
    geometry: MatrixGeometry,
) -> TractReadout {
    let MatrixGeometry {
        stride,
        read_len,
        reference_len,
        left_flank_len,
        right_flank_len,
        period,
        min_flank_support,
    } = geometry;
    let left_junction = left_flank_len;
    let right_junction = reference_len - right_flank_len;
    let mut tract_start = 0usize;
    let mut tract_end = read_len;
    let mut left_support = 0usize;
    let mut right_support = 0usize;

    let mut i = read_len;
    let mut j = reference_len;
    let mut state = final_state;
    while i != 0 || j != 0 {
        let pred = backpointers[i * stride + j][state.index()];
        match state {
            State::Match => {
                let consumed = j - 1;
                if consumed < left_junction {
                    left_support += 1;
                }
                if consumed >= right_junction {
                    right_support += 1;
                }
                if consumed == left_junction {
                    tract_start = i - 1;
                }
                if right_flank_len > 0 && consumed == right_junction {
                    tract_end = i - 1;
                }
                i -= 1;
                j -= 1;
            }
            State::Deletion => {
                let consumed = j - 1;
                if consumed == left_junction {
                    tract_start = i;
                }
                if right_flank_len > 0 && consumed == right_junction {
                    tract_end = i;
                }
                j -= 1;
            }
            State::Insertion => {
                i -= 1;
            }
            State::SlipInsertion => {
                i -= period;
            }
            State::SlipDeletion => {
                let deleted = (j - period)..j;
                if deleted.contains(&left_junction) {
                    tract_start = i;
                }
                if right_flank_len > 0 && deleted.contains(&right_junction) {
                    tract_end = i;
                }
                j -= period;
            }
        }
        state = pred;
    }

    // The support each flank must muster, clamped to the flank that exists — a repeat at a
    // contig edge is judged against its real flank rather than made unanchorable.
    let left_needed = min_flank_support.min(left_flank_len);
    let right_needed = min_flank_support.min(right_flank_len);

    let left_anchored = left_flank_len > 0 && tract_start != 0 && left_support >= left_needed;
    let right_anchored =
        right_flank_len > 0 && tract_end != read_len && right_support >= right_needed;

    // A flank that did not hold reports no boundary: the observed repeat runs to that end of
    // the read. Leaving the junction offset in place would make `length_lower_bound` report a
    // bound shorter than the bases actually seen.
    let tract_start = if left_anchored { tract_start } else { 0 };
    let tract_end = if right_anchored { tract_end } else { read_len };

    TractReadout {
        tract_start: tract_start as u64,
        tract_end: tract_end as u64,
        left_anchored,
        right_anchored,
    }
}

impl<E: Emission> BestPathAligner for SsrNoiseRobustAligner<E> {
    type Scratch = NoiseRobustScratch;
    type Output = RepeatSpan;
    type Context<'a> = RepeatContext<'a>;

    fn align(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: Self::Context<'_>,
        scratch: &mut Self::Scratch,
    ) -> Self::Output {
        self.delimit(read, reference, &context, scratch)
            .map_or(RepeatSpan::Unanchored, |readout| readout.classify())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::alignment::emission::PerQualityEmission;
    use crate::ng::alignment::{RepeatGeometry, StutterModel};
    use crate::ng::types::{Bp, Motif};

    /// The aperiodic flank bodies the synthetic bake-off uses, so the tests below sit in the
    /// same regime the configuration was chosen in.
    const BODY_L: &[u8] = b"GATCTTGCAAGCTGGAATCCGTTAC";
    const BODY_R: &[u8] = b"CAGTTCACGATCCTAAGGCTTGACG";

    fn frame(left: &[u8], tract: &[u8], right: &[u8], motif: &[u8]) -> (Vec<u8>, RepeatGeometry) {
        let mut reference = Vec::new();
        reference.extend_from_slice(left);
        reference.extend_from_slice(tract);
        reference.extend_from_slice(right);
        let geometry = RepeatGeometry {
            left_flank_len: Bp(left.len() as u64),
            right_flank_len: Bp(right.len() as u64),
            motif: Motif::new(motif).expect("a valid test motif"),
        };
        (reference, geometry)
    }

    fn repeat(motif: &[u8], units: usize) -> Vec<u8> {
        motif
            .iter()
            .copied()
            .cycle()
            .take(units * motif.len())
            .collect()
    }

    fn measure_with(
        config: NoiseRobustConfig,
        read: &[u8],
        reference: &[u8],
        geometry: &RepeatGeometry,
    ) -> RepeatSpan {
        let aligner = SsrNoiseRobustAligner::with_config(PerQualityEmission::new(), config);
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry,
            stutter: &stutter,
        };
        let quality = vec![40u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = NoiseRobustScratch::new();
        aligner.align(bases, reference, context, &mut scratch)
    }

    fn measure(read: &[u8], reference: &[u8], geometry: &RepeatGeometry) -> RepeatSpan {
        measure_with(NoiseRobustConfig::default(), read, reference, geometry)
    }

    /// The baseline sanity check: with no slip, the aligner must not invent one.
    #[test]
    fn a_clean_read_measures_the_reference_tract() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        assert_eq!(
            measure(&reference, &reference, &geometry).measured_length(),
            Some(24)
        );
    }

    /// Genuine length variation must still come through: the slip route is untouched by any
    /// of the three boundary changes.
    #[test]
    fn whole_unit_expansions_and_contractions_are_measured() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        for units in [5usize, 6, 7, 8, 9, 10, 12] {
            let mut read = BODY_L.to_vec();
            read.extend_from_slice(&repeat(b"AAG", units));
            read.extend_from_slice(BODY_R);
            assert_eq!(
                measure(&read, &reference, &geometry).measured_length(),
                Some(units as u64 * 3),
                "a {units}-unit allele was not measured at its own length"
            );
        }
    }

    /// **The out-of-frame route survives the junction guard.** A one-base insertion mid-tract
    /// is not a whole-unit slip, so it can only be spelled by the out-of-frame gap the guard
    /// narrows — and it must still be read out verbatim rather than rounded to a unit count.
    #[test]
    fn an_out_of_frame_insertion_mid_tract_is_still_measured() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        let mut tract = repeat(b"AAG", 8);
        tract.insert(12, b'T');
        let mut read = BODY_L.to_vec();
        read.extend_from_slice(&tract);
        read.extend_from_slice(BODY_R);
        assert_eq!(
            measure(&read, &reference, &geometry).measured_length(),
            Some(25)
        );
    }

    /// An interior substitution keeps the tract's own length — an interrupted repeat must not
    /// be pulled onto a tidy unit count.
    #[test]
    fn an_interior_substitution_keeps_the_tract_length() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        let mut tract = repeat(b"AAG", 8);
        tract[12] = b'T';
        let mut read = BODY_L.to_vec();
        read.extend_from_slice(&tract);
        read.extend_from_slice(BODY_R);
        assert_eq!(
            measure(&read, &reference, &geometry).measured_length(),
            Some(24)
        );
    }

    /// **The clamp is load-bearing.** On a short tract an unclamped six-base guard would close
    /// both ends onto each other and leave no column able to open an out-of-frame gap, which
    /// would force an impure short repeat onto a whole-unit length. Configuring an absurd guard
    /// must therefore change nothing: [`effective_guard`] cuts it back to what the tract can
    /// spare. The `junction_guard_bases: 0` run is the control — it shows the expected answer
    /// is reachable at all, so a pass here is the clamp working and not the assertion being
    /// vacuous.
    #[test]
    fn the_guard_is_clamped_so_a_short_tract_keeps_an_interior_route() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 4), BODY_R, b"AAG");
        let mut tract = repeat(b"AAG", 4);
        tract.insert(6, b'T');
        let mut read = BODY_L.to_vec();
        read.extend_from_slice(&tract);
        read.extend_from_slice(BODY_R);

        for junction_guard_bases in [0usize, 6, 100] {
            let config = NoiseRobustConfig {
                junction_guard_bases,
                ..NoiseRobustConfig::default()
            };
            assert_eq!(
                measure_with(config, &read, &reference, &geometry).measured_length(),
                Some(13),
                "a short impure tract lost its out-of-frame route at guard {junction_guard_bases}"
            );
        }
    }

    #[test]
    fn effective_guard_always_leaves_a_unit_of_tract_open() {
        // Long tract: the configured width stands.
        assert_eq!(effective_guard(6, 24, 3), 6);
        // Short tract: cut back to what is left after one unit, halved between the two ends.
        assert_eq!(effective_guard(6, 12, 3), 4);
        assert_eq!(effective_guard(6, 9, 3), 3);
        // Degenerate: a tract no longer than one unit can afford no guard at all, and the
        // arithmetic must saturate rather than wrap.
        assert_eq!(effective_guard(6, 3, 3), 0);
        assert_eq!(effective_guard(6, 0, 3), 0);
        assert_eq!(effective_guard(6, 1, 6), 0);
        // Every clamped width leaves at least one unit of columns openable.
        for tract_len in 0..40usize {
            for period in 1..=6usize {
                let guard = effective_guard(6, tract_len, period);
                assert!(
                    tract_len >= 2 * guard + period.min(tract_len),
                    "guard {guard} closed a tract of {tract_len} at period {period}"
                );
            }
        }
    }

    /// **A single stray base is not an anchor.** A read that stops one base into the right
    /// flank has crossed the junction, so algorithm 4's "did the traceback cross it" test
    /// calls the result a measurement — on the strength of one base, at a locus where the
    /// read may well have run off inside a longer allele. That is a long allele reported as a
    /// confident short one, the failure [`RepeatSpan`]'s partial variants exist to prevent.
    ///
    /// The `min_flank_support: 0` run reproduces that behaviour, so this shows the rule
    /// *changing* the answer rather than merely agreeing with it — and the sweep to five
    /// bases shows where the rule stops, since a five-base flank is a real anchor and must
    /// keep measuring.
    #[test]
    fn a_flank_held_by_one_base_does_not_anchor_a_measurement() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        let mut read = BODY_L.to_vec();
        read.extend_from_slice(&repeat(b"AAG", 8));
        read.extend_from_slice(&BODY_R[..1]);

        let unsupported = NoiseRobustConfig {
            min_flank_support: 0,
            ..NoiseRobustConfig::default()
        };
        assert_eq!(
            measure_with(unsupported, &read, &reference, &geometry).measured_length(),
            Some(24),
            "the case no longer reproduces the one-base anchor; the assertion below proves nothing"
        );

        let span = measure(&read, &reference, &geometry);
        assert!(
            matches!(span, RepeatSpan::FromLeft(_)),
            "a tract held by one flank base was still called a measurement: {span:?}"
        );
        // The bound must cover every base the read showed past the left junction, not stop at
        // the junction the discarded anchor claimed.
        assert_eq!(span.length_lower_bound(read.len() as u64), 25);
    }

    /// A read with plenty of flank on both sides still measures — the support rule must not
    /// turn ordinary reads into lower bounds. Five bases of flank in the read is production's
    /// own admission threshold, and the shipped support of three sits below it deliberately.
    #[test]
    fn a_short_but_real_flank_still_anchors() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        let mut read = BODY_L.to_vec();
        read.extend_from_slice(&repeat(b"AAG", 9));
        read.extend_from_slice(&BODY_R[..5]);
        assert_eq!(
            measure(&read, &reference, &geometry).measured_length(),
            Some(27)
        );
    }

    /// **Banning flank gaps must not close the partial routes.** A read that begins inside the
    /// tract crosses the left flank by deletion, and one that ends inside it crosses the right
    /// flank the same way; those crossings are exempt, and without the exemption every partial
    /// read would become a wrong answer rather than a lower bound.
    #[test]
    fn both_partial_directions_survive_the_flank_gap_ban() {
        let motif = b"AAG";
        let (reference, geometry) = frame(BODY_L, &repeat(motif, 8), BODY_R, motif);
        let long = repeat(motif, 48);

        // Runs off the 3′ end inside the tract.
        let mut from_left = BODY_L.to_vec();
        from_left.extend_from_slice(&long[..60]);
        assert!(matches!(
            measure(&from_left, &reference, &geometry),
            RepeatSpan::FromLeft(_)
        ));

        // Begins inside the tract and reaches the right flank.
        let mut from_right = long[long.len() - 60..].to_vec();
        from_right.extend_from_slice(BODY_R);
        assert!(matches!(
            measure(&from_right, &reference, &geometry),
            RepeatSpan::FromRight(_)
        ));
    }

    /// The degenerate input the trait contract names: an empty reference is not a locus, so the
    /// answer is `Unanchored` rather than an error or a panic.
    #[test]
    fn an_empty_reference_is_unanchored() {
        let geometry = RepeatGeometry {
            left_flank_len: Bp(0),
            right_flank_len: Bp(0),
            motif: Motif::new(b"AAG").expect("a valid motif"),
        };
        assert_eq!(measure(b"AAGAAG", b"", &geometry), RepeatSpan::Unanchored);
    }

    /// **Purity per call**, which the cohort's byte-identity guarantee rests on: a reused
    /// scratch must give the same answer as a fresh one, in any order.
    #[test]
    fn a_reused_scratch_changes_no_answer() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        let aligner = SsrNoiseRobustAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &stutter,
        };
        let mut reused = NoiseRobustScratch::new();
        // Longest first, so a stale row from a larger read would be visible in a shorter one.
        for units in [14usize, 3, 9, 1, 8] {
            let mut read = BODY_L.to_vec();
            read.extend_from_slice(&repeat(b"AAG", units));
            read.extend_from_slice(BODY_R);
            let quality = vec![40u8; read.len()];
            let bases = ReadBases::try_new(&read, &quality).expect("matched lengths");
            let fresh = aligner.align(bases, &reference, context, &mut NoiseRobustScratch::new());
            let again = aligner.align(bases, &reference, context, &mut reused);
            assert_eq!(
                fresh, again,
                "scratch reuse changed the answer at {units} units"
            );
        }
    }

    /// The shipped configuration, pinned. These three values are the finding the module
    /// records; a silent edit to any of them is a different algorithm wearing the same name.
    #[test]
    fn the_shipped_configuration_is_the_measured_one() {
        let config = NoiseRobustConfig::default();
        assert!(!config.flank_gaps);
        assert_eq!(config.junction_guard_bases, 6);
        assert_eq!(config.min_flank_support, 3);
        assert_eq!(config.slip_margin, 0.0);
        assert_eq!(
            SsrNoiseRobustAligner::new(PerQualityEmission::new()).config(),
            config
        );
    }
}
