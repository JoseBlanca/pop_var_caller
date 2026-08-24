//! Algorithm 4r — the **robust-indel** repeat delimiter: algorithm 4n with the flank gap
//! *priced* instead of *banned*.
//!
//! [`super::ssr_noise_robust`] (algorithm 4n) buys its noise resistance by closing the ordinary
//! per-base gap route outside the tract entirely. That works — a bad flank base is then priced
//! as exactly what it is, one mismatch, and no single base can slide the read against the frame
//! and move the junction. But a ban is not a price, and it says something false about biology: a
//! flank *can* carry a real 1 bp indel, and when it does, algorithm 4n has no way to spell it.
//! It must instead frameshift the rest of the flank against the reference, which mismeasures the
//! tract. On the synthetic bake-off that is the difference between `flank_indel = 0.000` and
//! `1.000`.
//!
//! **The fix is arithmetic, and the arithmetic is the whole design.** Two facts set the scale,
//! both at Q40:
//!
//! * A **real 1 bp flank indel** is enormously cheap to admit. The synthetic indel sits mid-flank,
//!   so refusing it forces roughly half the flank out of register — a dozen mismatches at
//!   `ln(1e-4/3) ≈ -10.31` each, on the order of **-120 nats**. One gap open buys all of that
//!   back. The gap route therefore has around a hundred nats of headroom before a genuine indel
//!   stops being worth taking.
//!
//! * A **noise-driven slide** is comparatively cheap to *deny*. Sliding the read one base against
//!   the frame and back needs an insertion **and** a deletion — two opens, `2 × ln(2.9e-5) ≈
//!   -20.9` — and it pays for itself as soon as it repairs three miscalled bases (`3 × 10.31 ≈
//!   30.9`). Three bad bases in a 25 bp flank is ordinary at the rates the bake-off sweeps, which
//!   is why algorithm 4 (unbanned, unpriced) loses a base of boundary under noise.
//!
//! The two live at very different scales, so a single number separates them.
//! [`RobustIndelConfig::flank_gap_open_penalty`] is charged on **each** flank open, so a slide —
//! which needs two — pays it twice while a genuine indel pays it once. That asymmetry is the
//! point: at a penalty of `P` nats a slide must now repair `(20.9 + 2P) / 10.31` bad bases to be
//! worth taking, so `P = 35` lifts the bar from three bad bases to nine, while a genuine indel
//! still clears its hundred-nat benefit with room to spare. **A flank gap becomes the explanation
//! of many bases rather than of one or two** — which is the rule algorithm 4n was reaching for,
//! stated as a price rather than as a prohibition.
//!
//! [`RobustIndelConfig::flank_junction_guard_bases`] closes the residual case the price alone
//! cannot reach. Within a base or two of the junction an indel is not a *flank* event at all — it
//! is indistinguishable from the tract being a base longer or shorter, and the tract's own
//! [`RobustIndelConfig::junction_guard_bases`] exists precisely so that reading is not made there.
//! Allowing a priced flank gap at the junction column would reopen from the outside the door the
//! tract guard shuts from the inside. So the innermost flank columns keep algorithm 4n's ban, and
//! the price governs everywhere else.
//!
//! Everything else is algorithm 4n's, unchanged: the five-state recurrence, the tract junction
//! guard, the minimum flank support, the terminal-route exemptions that keep
//! [`RepeatSpan::FromLeft`]/[`RepeatSpan::FromRight`] reachable, the tie-break order and the
//! traceback. In particular the **partial routes are untouched** — this algorithm only ever adds
//! paths a banned one lacked, and adding paths cannot turn a read that was a lower bound into no
//! answer at all.

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
/// aligner also uses it to *close* a route at the flank junction.
const UNREACHABLE: f64 = f64::NEG_INFINITY;

/// The five knobs that separate this delimiter from algorithm 4 — see the module docs.
///
/// It is **constructor state, not scratch**: each value changes results, so a bake-off has to be
/// able to see and report the setting it compared at.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RobustIndelConfig {
    /// Extra log-cost, in nats, charged on **every** ordinary gap opened outside the tract.
    /// `35.0` ships.
    ///
    /// This is algorithm 4n's `flank_gaps: bool` turned into a number, and the number is the
    /// module's finding. Zero reproduces algorithm 4's flank behaviour (a single noisy base can
    /// buy a slide); infinity reproduces algorithm 4n's ban (a real indel cannot be spelled).
    /// Thirty-five nats sits where a slide — which pays the penalty twice, once per open — needs
    /// roughly nine miscalled bases to be worth taking, while a genuine mid-flank indel, worth
    /// something like a dozen mismatches, still clears it paying once.
    ///
    /// The bake-off is flat across a wide band here (see the sweep recorded on
    /// [`Self::default`]), which is the reassuring shape: the two regimes really are an order of
    /// magnitude apart, so the setting is a plateau rather than a knife edge.
    pub flank_gap_open_penalty: f64,
    /// How many flank columns adjacent to each junction may not open an ordinary gap at all.
    /// `2` ships.
    ///
    /// Priced or not, a gap in the first flank column is not a flank event: it is the tract being
    /// one base longer or shorter, spelled from the outside. [`Self::junction_guard_bases`] shuts
    /// that reading from the tract side, and this shuts it from the flank side; leaving one open
    /// would make the other pointless. Kept deliberately narrow — two bases, not the tract
    /// guard's six — because a real indel one base past the junction is still a real indel, and
    /// the price, not the ban, is meant to be doing the work everywhere it can.
    pub flank_junction_guard_bases: usize,
    /// How many **tract** columns at each end are barred from opening an out-of-frame gap.
    /// `6` ships — algorithm 4n's value and rationale, unchanged.
    ///
    /// Always clamped so the tract keeps an interior: see [`effective_guard`].
    pub junction_guard_bases: usize,
    /// Extra log-cost per whole-unit slip run, in nats. `0.0` ships — algorithm 4n's finding,
    /// unchanged: the slip route is already priced by the stutter model, and taxing it further
    /// costs more on genuine expansions than it recovers on noise.
    pub slip_margin: f64,
    /// How many read bases must be **aligned against a flank** before that flank counts as
    /// anchoring the tract. `3` ships — algorithm 4n's value and rationale, unchanged.
    ///
    /// Counted in aligned columns, not in agreeing bases, so a noisy flank still anchors while a
    /// flank represented by one lone base does not. Clamped to the flank that actually exists.
    pub min_flank_support: usize,
}

impl Default for RobustIndelConfig {
    /// The shipped configuration, chosen on the synthetic bake-off.
    ///
    /// The penalty sweep is the reason this algorithm exists, so it is worth recording rather
    /// than merely asserting. Holding everything else fixed and varying
    /// [`Self::flank_gap_open_penalty`] moved the composite as:
    ///
    /// ```text
    /// penalty   clean  partial  p_noise    noise  flank_indel  COMPOSITE
    ///    0.0    1.000    1.000  0.95637  0.99657        1.000   0.988235
    ///    5.0    1.000    1.000  0.97673  0.99883        1.000   0.993890
    ///   10.0    1.000    1.000  0.98460  0.99895        1.000   0.995888
    ///   15.0    1.000    1.000  0.98637  0.99907        1.000   0.996360
    ///   20.0    1.000    1.000  0.98687  0.99907        1.000   0.996487
    ///   25.0    1.000    1.000  0.98705  0.99909        1.000   0.996536   <- ceiling reached
    ///   35.0    1.000    1.000  0.98705  0.99909        1.000   0.996536   <- shipped
    ///   50.0    1.000    1.000  0.98705  0.99909        1.000   0.996536
    ///   55.0    1.000    1.000  0.98705  0.99909        0.833   0.954869   <- cliff
    ///   90.0    1.000    1.000  0.98705  0.99909        0.375   0.840300
    ///    inf    1.000    1.000  0.98705  0.99909        0.000   0.746536   (= algorithm 4n)
    /// ```
    ///
    /// `noise` and `p_noise` climb to algorithm 4n's own ceiling by 25 nats and stop there;
    /// `flank_indel` holds at 1.000 until about 52, where the gap finally becomes dearer than the
    /// frameshift it was bought to avoid, and then falls away to the ban. The usable band is
    /// therefore `[25, 50]`, and 35 is taken from the middle of it. The two edges are not alike —
    /// the lower one degrades a thousandth at a time while the upper one is a cliff — so the
    /// choice leans slightly low, which is the side that fails gracefully.
    fn default() -> Self {
        Self {
            flank_gap_open_penalty: 35.0,
            flank_junction_guard_bases: 2,
            junction_guard_bases: 6,
            slip_margin: 0.0,
            min_flank_support: 3,
        }
    }
}

/// How wide the **tract** junction guard may actually be at this locus: the configured width, cut
/// back so the tract keeps at least one motif unit of columns that can still open an out-of-frame
/// gap. Algorithm 4n's clamp, unchanged — without it a short tract would lose its whole interior
/// route and every interrupted repeat there would be forced onto a tidy unit count.
#[inline]
fn effective_guard(configured: usize, tract_len: usize, period: usize) -> usize {
    configured.min(tract_len.saturating_sub(period) / 2)
}

/// The whole-unit slip transition costs, in log space, derived from the shared [`StutterModel`].
#[derive(Debug, Clone, Copy)]
struct SlipCosts {
    /// `ln(whole_repeat_longer_share · whole_repeat_one_step_share)
    /// − ln(same_length_share) − margin`.
    open_expansion: f64,
    /// `ln(whole_repeat_shorter_share · whole_repeat_one_step_share)
    /// − ln(same_length_share) − margin`.
    open_contraction: f64,
    /// `ln(1 − whole_repeat_one_step_share)`, per unit after the first.
    extend: f64,
}

impl SlipCosts {
    fn from_model(model: &StutterModel, margin: f64) -> Self {
        let ln_same_length_share = model.same_length_share().ln();
        Self {
            open_expansion: (model.whole_repeat_longer_share()
                * model.whole_repeat_one_step_share())
            .ln()
                - ln_same_length_share
                - margin,
            open_contraction: (model.whole_repeat_shorter_share()
                * model.whole_repeat_one_step_share())
            .ln()
                - ln_same_length_share
                - margin,
            extend: (1.0 - model.whole_repeat_one_step_share()).ln(),
        }
    }
}

/// Per-worker scratch: a ring of `period + 1` score rows plus the full backpointer matrix, as the
/// whole-unit states need (a unit insertion reaches `(i, j)` from `(i − period, j)`).
///
/// Grow-and-keep; buffers only, deciding nothing that changes a result.
#[derive(Debug, Default)]
pub struct RobustIndelScratch {
    rows: Vec<Vec<[f64; STATES]>>,
    backpointers: Vec<[State; STATES]>,
}

impl RobustIndelScratch {
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

/// Pick the best candidate, keeping the **first on ties** — the spec's determinism rule, encoded
/// by the caller's argument order (match, then deletion, then insertion, then slips).
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

/// Algorithm 4r: the robust-indel repeat delimiter. Generic over its [`Emission`], never behind
/// `dyn` (arch §4).
#[derive(Debug, Clone, Copy)]
pub struct SsrRobustIndelAligner<E: Emission> {
    emission: E,
    costs: TransitionCosts,
    config: RobustIndelConfig,
}

impl<E: Emission> SsrRobustIndelAligner<E> {
    /// Build the delimiter at its shipped configuration.
    #[must_use]
    pub fn new(emission: E) -> Self {
        Self::with_config(emission, RobustIndelConfig::default())
    }

    /// Build the delimiter at an explicit configuration — the form a bake-off uses, so the setting
    /// under comparison is named at the call site.
    #[must_use]
    pub fn with_config(emission: E, config: RobustIndelConfig) -> Self {
        Self {
            emission,
            costs: TransitionCosts::new(),
            config,
        }
    }

    /// The configuration this delimiter is running at.
    #[must_use]
    pub fn config(&self) -> RobustIndelConfig {
        self.config
    }

    /// Measure the read's repeat. `None` when there is no reference frame — a defined answer, not
    /// an error (arch §3).
    #[must_use]
    #[allow(clippy::too_many_lines)]
    pub fn delimit(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: &RepeatContext<'_>,
        scratch: &mut RobustIndelScratch,
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

        let guard = effective_guard(
            self.config.junction_guard_bases,
            tract_last - left_flank_len,
            period,
        );
        let flank_guard = self.config.flank_junction_guard_bases;
        // The **priced** flank open: the ordinary open, made dearer by the configured penalty.
        // Charged per open, so a slide (an insertion plus a deletion) pays it twice.
        let flank_gap_open = self.costs.ln_gap_open() - self.config.flank_gap_open_penalty;

        // The per-base gap-open cost by column: the tract-grade open inside the tract,
        // `UNREACHABLE` at the tract's own junction guard, `UNREACHABLE` in the innermost flank
        // columns, and the priced open everywhere else in the flank.
        let gap_open = |column: usize| {
            if column_in_tract(column) {
                let near_junction = column <= left_flank_len + guard || column + guard > tract_last;
                if guard > 0 && near_junction {
                    UNREACHABLE
                } else {
                    self.costs.ln_gap_open_tract()
                }
            } else {
                // Distance, in flank columns, from the junction this column sits against.
                let distance = if column <= left_flank_len {
                    left_flank_len - column
                } else {
                    column - tract_last - 1
                };
                if distance < flank_guard {
                    UNREACHABLE
                } else {
                    flank_gap_open
                }
            }
        };

        // The **terminal** deletion route, exempt from every flank rule: a read that stops inside
        // the tract still has to reach the frame's end, and one that starts inside it still has to
        // be reached from the frame's start. Those are the `FromLeft`/`FromRight` cases, and
        // pricing them would turn every partial read into a wrong answer.
        let terminal_gap_open = self.costs.ln_gap_open();

        // The `k`-th base of a whole unit inserted at tract column `column`: the motif,
        // phase-aligned so a unit inserted at a unit boundary begins at `motif[0]`.
        let motif_base = |column: usize, k: usize| motif[((column - left_flank_len) + k) % period];

        let RobustIndelScratch { rows, backpointers } = scratch;

        // Row 0 — no read base consumed. Every deletion here is a *leading* one (the read starts
        // later in the frame), so it takes the terminal cost regardless of the flank rules.
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
            // boundary-sliding gap the flank price is aimed at.
            let trailing = row_index == read_len;

            // Column 0 — a read base before any reference base: the leading-insertion route, which
            // is how a read that overhangs the frame's start is placed. Kept at the ordinary open
            // cost for the same reason the terminal deletions are.
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
                // last read row the terminal route applies, so a priced or closed flank open
                // reverts to the ordinary one.
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
    /// [`RobustIndelConfig::min_flank_support`], carried into the walk that applies it.
    min_flank_support: usize,
}

/// Walk the traceback, read the tract off the two flank junctions, and check that each flank
/// actually held it — algorithm 4n's walk, unchanged.
///
/// The support counts are read bases *aligned against* a flank column, not bases that agree with
/// it: a miscalled flank base still holds the anchor, while a flank represented by one lone read
/// base does not.
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

    // The support each flank must muster, clamped to the flank that exists — a repeat at a contig
    // edge is judged against its real flank rather than made unanchorable.
    let left_needed = min_flank_support.min(left_flank_len);
    let right_needed = min_flank_support.min(right_flank_len);

    let left_anchored = left_flank_len > 0 && tract_start != 0 && left_support >= left_needed;
    let right_anchored =
        right_flank_len > 0 && tract_end != read_len && right_support >= right_needed;

    // A flank that did not hold reports no boundary: the observed repeat runs to that end of the
    // read. Leaving the junction offset in place would make `length_lower_bound` report a bound
    // shorter than the bases actually seen.
    let tract_start = if left_anchored { tract_start } else { 0 };
    let tract_end = if right_anchored { tract_end } else { read_len };

    TractReadout {
        tract_start: tract_start as u64,
        tract_end: tract_end as u64,
        left_anchored,
        right_anchored,
    }
}

impl<E: Emission> BestPathAligner for SsrRobustIndelAligner<E> {
    type Scratch = RobustIndelScratch;
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
    use crate::ng::alignment::{RepeatGeometry, StutterModel, StutterRates};
    use crate::ng::types::{Bp, Motif};

    /// The aperiodic flank bodies the synthetic bake-off uses, so the tests below sit in the same
    /// regime the configuration was chosen in.
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
        config: RobustIndelConfig,
        read: &[u8],
        reference: &[u8],
        geometry: &RepeatGeometry,
    ) -> RepeatSpan {
        let aligner = SsrRobustIndelAligner::with_config(PerQualityEmission::new(), config);
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry,
            stutter: &stutter,
        };
        let quality = vec![40u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = RobustIndelScratch::new();
        aligner.align(bases, reference, context, &mut scratch)
    }

    fn measure(read: &[u8], reference: &[u8], geometry: &RepeatGeometry) -> RepeatSpan {
        measure_with(RobustIndelConfig::default(), read, reference, geometry)
    }

    /// Contraction-biased parameters, with **all six rates different** — HipSTR's fitted
    /// values are contraction-biased, and a fixture that gives the two members of a pair the
    /// same value cannot tell them apart.
    ///
    /// This module's alignment fixtures use `StutterModel::hipstr_shipped()`, whose two
    /// whole-repeat shares are both 0.05 — so `open_expansion` and `open_contraction` derive
    /// to the *same number*, and swapping them was a change no test in this file could see.
    /// Measured: both opened at −2.9191921964316565 clean and swapped, bit for bit. The
    /// cost-level tests below use this fixture instead, which is the same discipline
    /// `all_distinct()` enforces in `stutter.rs`.
    fn contraction_biased() -> StutterModel {
        StutterModel::new(StutterRates {
            whole_repeat_longer_share: 0.03,
            whole_repeat_shorter_share: 0.07,
            whole_repeat_one_step_share: 0.9,
            part_repeat_longer_share: 0.004,
            part_repeat_shorter_share: 0.012,
            part_repeat_one_step_share: 0.8,
        })
    }

    /// **A contraction must be cheaper to open than an expansion**, on a model that says so.
    /// The asymmetry lives at the cost level, and this is the only thing in this file that
    /// reads it: without it, `open_expansion` and `open_contraction` could be exchanged and
    /// every other test here would still pass.
    #[test]
    fn a_contraction_is_cheaper_to_open_than_an_expansion() {
        let slip = SlipCosts::from_model(&contraction_biased(), 0.0);
        assert!(
            slip.open_contraction > slip.open_expansion,
            "a contraction-biased model must make contraction the cheaper slip to open"
        );
    }

    /// The slip costs reconstruct the stutter model's own probabilities exactly — the affine
    /// open plus its extends must sum back to `ln(P(n)) − ln(same_length_share)`, or this
    /// aligner is pricing a different distribution than the one it was handed. Run at a zero margin, since a
    /// non-zero one shifts both opens by the same constant and so cannot be checked against
    /// the model's own probabilities; `a_contraction_is_cheaper_to_open_than_an_expansion`
    /// covers the pairing at any margin.
    #[test]
    fn the_slip_costs_reconstruct_the_stutter_probability() {
        let model = contraction_biased();
        let slip = SlipCosts::from_model(&model, 0.0);
        let period = std::num::NonZeroU8::new(3).unwrap();
        let ln_same_length_share = model.same_length_share().ln();

        for n in 1..=5i64 {
            // Expansion of n units: one open plus (n − 1) extends, relative to no slip.
            let reconstructed = slip.open_expansion + (n - 1) as f64 * slip.extend;
            let expected = model.probability(n * 3, period).ln() - ln_same_length_share;
            assert!(
                (reconstructed - expected).abs() < 1e-12,
                "expansion of {n} units diverged: {reconstructed} vs {expected}"
            );
            // Contraction likewise, and it is the direction the fixture makes cheaper.
            let reconstructed = slip.open_contraction + (n - 1) as f64 * slip.extend;
            let expected = model.probability(-n * 3, period).ln() - ln_same_length_share;
            assert!(
                (reconstructed - expected).abs() < 1e-12,
                "contraction of {n} units diverged"
            );
        }
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

    /// Genuine length variation must still come through: the slip route is untouched.
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

    /// **This is the algorithm's reason for existing.** A real 1 bp indel in the read's flank must
    /// not disturb the tract measurement — and it can only be spelled by the very gap algorithm 4n
    /// bans. Both flanks, insertion and deletion.
    ///
    /// The `f64::INFINITY` run reproduces algorithm 4n's ban, so this shows the price *changing*
    /// the answer rather than merely agreeing with it.
    #[test]
    fn a_one_base_flank_indel_leaves_the_tract_measurement_alone() {
        let motif = b"AAG";
        let (reference, geometry) = frame(BODY_L, &repeat(motif, 8), BODY_R, motif);
        let banned = RobustIndelConfig {
            flank_gap_open_penalty: f64::INFINITY,
            ..RobustIndelConfig::default()
        };

        for right_side in [false, true] {
            for insert in [false, true] {
                let edit = |flank: &[u8]| -> Vec<u8> {
                    let mut f = flank.to_vec();
                    let at = f.len() / 2;
                    if insert {
                        f.insert(at, b'T');
                    } else {
                        f.remove(at);
                    }
                    f
                };
                let (left, right) = if right_side {
                    (BODY_L.to_vec(), edit(BODY_R))
                } else {
                    (edit(BODY_L), BODY_R.to_vec())
                };
                let mut read = left;
                read.extend_from_slice(&repeat(motif, 8));
                read.extend_from_slice(&right);

                assert_eq!(
                    measure(&read, &reference, &geometry).measured_length(),
                    Some(24),
                    "a 1 bp flank indel (right_side={right_side}, insert={insert}) moved the tract"
                );
                assert_ne!(
                    measure_with(banned, &read, &reference, &geometry).measured_length(),
                    Some(24),
                    "the banned control no longer fails, so the assertion above proves nothing \
                     (right_side={right_side}, insert={insert})"
                );
            }
        }
    }

    /// **The price must not be a free pass.** A flank gap has to stay dearer than the handful of
    /// mismatches it would otherwise repair, which is what keeps a noisy base from sliding the
    /// boundary. Pricing it is only a defensible middle ground if the shipped penalty really is
    /// worth several mismatches at Q40 — otherwise this is algorithm 4 again.
    #[test]
    fn the_shipped_penalty_is_worth_more_than_two_mismatches() {
        let mismatch = (1e-4f64 / 3.0).ln().abs();
        let penalty = RobustIndelConfig::default().flank_gap_open_penalty;
        assert!(
            penalty > 2.0 * mismatch,
            "a flank gap priced at {penalty} nats is cheaper than the two mismatches it must beat"
        );
    }

    /// The out-of-frame route survives the tract junction guard: a one-base insertion mid-tract is
    /// not a whole-unit slip, so it can only be spelled by the gap the guard narrows — and it must
    /// still be read out verbatim rather than rounded to a unit count.
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

    /// An interior substitution keeps the tract's own length — an interrupted repeat must not be
    /// pulled onto a tidy unit count.
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

    /// The tract guard clamp is load-bearing: on a short tract an unclamped six-base guard would
    /// leave no column able to open an out-of-frame gap, forcing an impure short repeat onto a
    /// whole-unit length. The `junction_guard_bases: 0` run is the control.
    #[test]
    fn the_guard_is_clamped_so_a_short_tract_keeps_an_interior_route() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 4), BODY_R, b"AAG");
        let mut tract = repeat(b"AAG", 4);
        tract.insert(6, b'T');
        let mut read = BODY_L.to_vec();
        read.extend_from_slice(&tract);
        read.extend_from_slice(BODY_R);

        for junction_guard_bases in [0usize, 6, 100] {
            let config = RobustIndelConfig {
                junction_guard_bases,
                ..RobustIndelConfig::default()
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
        assert_eq!(effective_guard(6, 24, 3), 6);
        assert_eq!(effective_guard(6, 12, 3), 4);
        assert_eq!(effective_guard(6, 9, 3), 3);
        assert_eq!(effective_guard(6, 3, 3), 0);
        assert_eq!(effective_guard(6, 0, 3), 0);
        assert_eq!(effective_guard(6, 1, 6), 0);
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

    /// A single stray base is not an anchor: a read that stops one base into the right flank must
    /// come back a lower bound, not a confident short measurement. The `min_flank_support: 0` run
    /// reproduces the failure, so the rule is shown changing the answer.
    #[test]
    fn a_flank_held_by_one_base_does_not_anchor_a_measurement() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        let mut read = BODY_L.to_vec();
        read.extend_from_slice(&repeat(b"AAG", 8));
        read.extend_from_slice(&BODY_R[..1]);

        let unsupported = RobustIndelConfig {
            min_flank_support: 0,
            ..RobustIndelConfig::default()
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
        assert_eq!(span.length_lower_bound(read.len() as u64), 25);
    }

    /// A read with plenty of flank on both sides still measures — the support rule must not turn
    /// ordinary reads into lower bounds.
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

    /// **Pricing the flank gap must not cost a partial its lower bound.** A read that begins inside
    /// the tract crosses the left flank by deletion, and one that ends inside it crosses the right
    /// flank the same way; those crossings take the terminal cost, never the priced one. Turning a
    /// partial into `Unanchored` destroys real evidence, so this is the invariant that outranks
    /// every score on the bake-off.
    #[test]
    fn both_partial_directions_survive_the_flank_gap_price() {
        let motif = b"AAG";
        let (reference, geometry) = frame(BODY_L, &repeat(motif, 8), BODY_R, motif);
        let long = repeat(motif, 48);

        let mut from_left = BODY_L.to_vec();
        from_left.extend_from_slice(&long[..60]);
        assert!(matches!(
            measure(&from_left, &reference, &geometry),
            RepeatSpan::FromLeft(_)
        ));

        let mut from_right = long[long.len() - 60..].to_vec();
        from_right.extend_from_slice(BODY_R);
        assert!(matches!(
            measure(&from_right, &reference, &geometry),
            RepeatSpan::FromRight(_)
        ));
    }

    /// **The invariant that outranks the bake-off**, stated as a sweep rather than as a single
    /// case: under heavy miscalling a read that runs off inside a long tract must come back a
    /// lower bound. It may not come back `Unanchored`.
    ///
    /// The direction matters and only one way is allowed. Turning an over-claimed *complete* into
    /// a partial is this family's whole purpose — it replaces a wrong number with an honest one.
    /// Turning a *partial* into `Unanchored` throws the read away, and a downstream genotype
    /// cannot recover evidence that was never reported. The flank gap price is a new route into
    /// the flank, so it is exactly the kind of change that could cost a partial its anchor, and it
    /// is checked at four error rates well past anything a real run sees.
    #[test]
    fn a_noisy_partial_never_degrades_to_no_answer() {
        let motif = b"AAG";
        let (reference, geometry) = frame(BODY_L, &repeat(motif, 8), BODY_R, motif);
        let long = repeat(motif, 48);
        // A deterministic bit-mixer, so the sweep is reproducible without a dependency.
        let mut state = 0x2545_F491_4F6C_DD1Du64;
        let mut next = move || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state
        };

        let mut from_left = BODY_L.to_vec();
        from_left.extend_from_slice(&long[..60]);
        let mut from_right = long[long.len() - 60..].to_vec();
        from_right.extend_from_slice(BODY_R);

        for template in [&from_left, &from_right] {
            for miscall_in in [200u64, 100, 50, 25] {
                for _ in 0..200 {
                    let read: Vec<u8> = template
                        .iter()
                        .map(|&b| {
                            if next() % miscall_in == 0 {
                                *b"ACGT".iter().find(|&&c| c != b).expect("a different base")
                            } else {
                                b
                            }
                        })
                        .collect();
                    let span = measure(&read, &reference, &geometry);
                    assert!(
                        !matches!(span, RepeatSpan::Unanchored),
                        "a noisy partial was thrown away at a 1-in-{miscall_in} miscall rate"
                    );
                }
            }
        }
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

    /// **Purity per call**, which the cohort's byte-identity guarantee rests on: a reused scratch
    /// must give the same answer as a fresh one, in any order.
    #[test]
    fn a_reused_scratch_changes_no_answer() {
        let (reference, geometry) = frame(BODY_L, &repeat(b"AAG", 8), BODY_R, b"AAG");
        let aligner = SsrRobustIndelAligner::new(PerQualityEmission::new());
        let stutter = StutterModel::hipstr_shipped();
        let context = RepeatContext {
            geometry: &geometry,
            stutter: &stutter,
        };
        let mut reused = RobustIndelScratch::new();
        for units in [14usize, 3, 9, 1, 8] {
            let mut read = BODY_L.to_vec();
            read.extend_from_slice(&repeat(b"AAG", units));
            read.extend_from_slice(BODY_R);
            let quality = vec![40u8; read.len()];
            let bases = ReadBases::try_new(&read, &quality).expect("matched lengths");
            let fresh = aligner.align(bases, &reference, context, &mut RobustIndelScratch::new());
            let again = aligner.align(bases, &reference, context, &mut reused);
            assert_eq!(
                fresh, again,
                "scratch reuse changed the answer at {units} units"
            );
        }
    }

    /// The shipped configuration, pinned. These values are the finding the module records; a
    /// silent edit to any of them is a different algorithm wearing the same name.
    #[test]
    fn the_shipped_configuration_is_the_measured_one() {
        let config = RobustIndelConfig::default();
        assert_eq!(config.flank_gap_open_penalty, 35.0);
        assert_eq!(config.flank_junction_guard_bases, 2);
        assert_eq!(config.junction_guard_bases, 6);
        assert_eq!(config.min_flank_support, 3);
        assert_eq!(config.slip_margin, 0.0);
        assert_eq!(
            SsrRobustIndelAligner::new(PerQualityEmission::new()).config(),
            config
        );
    }
}
