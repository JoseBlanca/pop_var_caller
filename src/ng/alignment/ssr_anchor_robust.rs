//! Algorithm 5r — **anchor-robust**: algorithm 5's firm anchor test, with the boundary itself
//! made resistant to single sequencing errors, and *without* algorithm 4n's ban on flank gaps.
//!
//! Three delimiters meet here, and it is worth being precise about what each contributes,
//! because the composite score is the sum of three separate failures being fixed:
//!
//! * **Algorithm 4** ([`super::ssr_best_path_unit_slip`]) supplies the recurrence: a five-state
//!   Viterbi whose whole-unit slips are priced from the shared [`StutterModel`].
//! * **Algorithm 5** ([`super::ssr_anchor_firm`]) supplies the *classification*: a flank counts
//!   as anchored only when the read actually **matched its bases** near the junction, which is
//!   what stops a long allele that outruns the read from being reported as a fabricated
//!   complete length.
//! * **Algorithm 4n** ([`super::ssr_noise_robust`]) supplies the idea this algorithm borrows
//!   for its noise resistance: near a junction, a single miscalled base must not be able to
//!   buy a one-base slide of the read against the reference.
//!
//! # The arithmetic the noise fix answers
//!
//! At Q40 a mismatch scores `ln(1e-4/3) ≈ -10.31`; the flank gap-open scores
//! `ln(2.9e-5) ≈ -10.45`. **One bad base costs about as much as opening a gap.** An insertion
//! plus a deletion slide the read one base against the reference for roughly the price of two
//! mismatches, so any noisy read where a one-base slide repairs three miscalls prefers the
//! slide — and the slide moves the tract junction. The measurement then comes back wrong by a
//! base, or by a whole unit once the slide lets a slip re-phase. That is the entire gap
//! between algorithm 5's `noise` axis and algorithm 4n's.
//!
//! # What this algorithm changes, and what it deliberately does not
//!
//! **[`JunctionGuard::tract_bases`] — the borrowed mechanism.** The out-of-frame per-base gap
//! may not *open* in the first or last few tract columns. Away from the junction that route is
//! doing its real job (spelling a genuinely impure repeat); beside the junction it is doing
//! something else entirely — trading one mismatch for a one-base boundary shift. Whole-unit
//! slips are untouched and may still open anywhere in the tract, so real length variation is
//! unaffected; only the sub-unit nudge at the boundary goes. Gap *extends* are untouched too,
//! which is what keeps a partial read able to cross the rest of the frame.
//!
//! **[`JunctionGuard::flank_bases`] — the same idea, from the other side.** A slide needs a gap
//! on *either* side of the junction, so guarding only the tract half leaves the flank half
//! open. Algorithm 4n closes it by banning flank gaps outright, which costs it the entire
//! flank-indel axis: a read carrying a real 1 bp indel in its flank then has no way to spell
//! it, and the tract measurement goes with it. This algorithm bans the gap-*open* only in the
//! handful of flank columns abutting the junction — the ones that can move a boundary — and
//! leaves the rest of the flank free. A real flank indel sits in the body of the flank, where
//! it is still representable; a noise-driven slide sits at the junction, where it is not.
//!
//! **Both guards are clamped** ([`JunctionGuard::resolve`]) so the tract keeps a motif unit of
//! openable interior and each flank keeps at least half of itself openable. The clamp is the
//! difference between a guard and a ban: a short tract must not have its whole interior closed,
//! or an interrupted repeat there is forced onto a tidy unit count — the exact failure the
//! out-of-frame route exists to prevent.
//!
//! **The terminal routes are exempt.** A read that stops inside the tract still has to reach
//! the frame's end, and one that starts inside it still has to be reached from the frame's
//! start; both cross guarded columns by deletion. Those crossings are what
//! [`RepeatSpan::FromLeft`] / [`RepeatSpan::FromRight`] report, so a deletion on read row 0
//! (nothing consumed yet) or on the last read row (nothing left to consume) keeps the ordinary
//! open cost. They are recognised structurally, never by a threshold.
//!
//! # The classification: demote a measurement, never a bound
//!
//! Algorithm 5's anchor quorum is kept in shape — [`AnchorRule`], the reference bases abutting
//! each junction, some of which must have *agreed* for the flank to count — and changed in two
//! ways.
//!
//! **The quorum is tightened**, from two agreements in five bases to four in seven, and this is
//! the single largest term in the score. The reasoning is in [`AnchorRule::default`]: a
//! fabricated crossing agrees with the flank only by chance, so raising the bar catches five
//! sixths of the ones algorithm 5 waves through, while a genuine read at a realistic miscall
//! rate clears it essentially always.
//!
//! **The floor is added**, and it is what makes the tightening safe. Algorithm 5 applies the
//! quorum to each side independently, so
//! a read whose path crossed both junctions but met the quorum on neither falls all the way
//! through to [`RepeatSpan::Unanchored`] — no observation at all. That is a demotion too far.
//! A read that reached a junction has *located* itself in the frame even when its flank bases
//! are too noisy to certify a length, and a lower bound is a usable observation where none is
//! not.
//!
//! So the rule here is stated as a floor rather than as two independent tests
//! (see [`Anchoring::decide`]):
//!
//! 1. The quorum may demote an over-claimed **complete** read to a lower bound.
//! 2. It may **never** demote a lower bound to nothing. A read whose path crossed only one
//!    junction keeps that junction, whatever the quorum says.
//! 3. A complete read that fails the quorum on *both* sides keeps the better-evidenced side as
//!    a bound, rather than being discarded.
//!
//! Under noise this only ever converts a wrong answer into a weaker one, and on real data it is
//! the difference between halving the partial yield and not.

use super::emission::Emission;
use super::stutter::StutterModel;
use super::{BestPathAligner, ReadBases, RepeatContext, RepeatSpan};

// Algorithm 3's tract-aware per-base gap transitions, shared rather than re-derived.
use super::ssr_best_path_flat_gap::TransitionCosts;

/// The five states a cell can be entered in — algorithm 4's, unchanged.
///
/// `#[repr(u8)]` with explicit discriminants so the values are stable array indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
enum State {
    /// A read base and a reference base consumed together.
    Match = 0,
    /// One read base against no reference base — an *out-of-frame* insertion.
    Insertion = 1,
    /// One reference base against no read base — an *out-of-frame* deletion.
    Deletion = 2,
    /// A whole unit of read bases the reference does not have — an *in-frame* expansion.
    SlipInsertion = 3,
    /// A whole unit of reference bases the read does not have — an *in-frame* contraction.
    SlipDeletion = 4,
}

/// The number of states — the width of each cell's score array and backpointer array.
const STATES: usize = 5;

impl State {
    #[inline]
    const fn index(self) -> usize {
        self as usize
    }
}

/// An unreachable score. A real value: some cells cannot be entered in some states, and this
/// aligner also uses it to *close* a route — a guarded gap-open is an unreachable transition,
/// not merely an expensive one.
const UNREACHABLE: f64 = f64::NEG_INFINITY;

/// How wide the no-gap-open guard is on each side of each junction.
///
/// **Constructor state, not scratch**: both widths change results, so an experiment has to be
/// able to see and report the setting it compared at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JunctionGuard {
    /// How many *tract* columns at each end may not open an out-of-frame gap.
    pub tract_bases: usize,
    /// How many *flank* columns abutting the tract may not open an out-of-frame gap. Keeping
    /// this well short of the flank length is what preserves a real flank indel's route.
    pub flank_bases: usize,
}

impl Default for JunctionGuard {
    /// **One motif unit on each side of each junction** — six bases, the widest period a
    /// [`super::Motif`] can hold, so the guard covers a whole unit's worth of re-phasing at
    /// every period and the rule is the same on both sides.
    ///
    /// The bake-off is flat across a wide neighbourhood here (anything from three to twelve
    /// tract columns, and four to twenty-five flank columns, scores within one part in ten
    /// thousand), so the width is chosen for the reason it can be stated, not for a decimal.
    /// What the flank width must *not* do is reach the body of the flank, where real indels
    /// live — six leaves nineteen of a twenty-five base flank openable.
    fn default() -> Self {
        Self {
            tract_bases: 6,
            flank_bases: 6,
        }
    }
}

impl JunctionGuard {
    /// The widths this locus can actually afford, given its tract and flanks.
    ///
    /// The tract keeps at least one motif unit of openable interior; each flank keeps at least
    /// half of itself openable. Without the clamp a short tract or a stub flank would have its
    /// whole out-of-frame route closed, which is a ban, not a guard.
    #[must_use]
    fn resolve(self, tract_len: usize, left_flank_len: usize, right_flank_len: usize) -> Resolved {
        Resolved {
            tract: self
                .tract_bases
                .min(tract_len.saturating_sub(WIDEST_PERIOD) / 2),
            left_flank: self.flank_bases.min(left_flank_len / 2),
            right_flank: self.flank_bases.min(right_flank_len / 2),
        }
    }
}

/// The interior a tract must keep openable, in bases. A whole motif unit is what it should be,
/// but the period is not in scope at [`JunctionGuard::resolve`]; the widest period a
/// [`super::Motif`] can hold is the conservative reading of "a unit", and on any tract long
/// enough for the guard to bind it is the same answer.
const WIDEST_PERIOD: usize = 6;

/// The per-locus guard widths, after clamping.
#[derive(Debug, Clone, Copy)]
struct Resolved {
    tract: usize,
    left_flank: usize,
    right_flank: usize,
}

/// The whole-unit slip transition costs, in log space, derived **from the shared
/// [`StutterModel`]** — not a second copy of its parameters.
///
/// `open` carries the affine open plus the `− ln(same_length_share)` baseline shift (a
/// best-path aligner maximises, so only scores relative to "no slip" matter); `extend` is the
/// geometric's per-extra-unit factor.
#[derive(Debug, Clone, Copy)]
struct SlipCosts {
    /// `ln(whole_repeat_longer_share · whole_repeat_one_step_share) − ln(same_length_share)`.
    open_expansion: f64,
    /// `ln(whole_repeat_shorter_share · whole_repeat_one_step_share) − ln(same_length_share)`.
    open_contraction: f64,
    /// `ln(1 − whole_repeat_one_step_share)`, charged per unit after the first.
    extend: f64,
}

impl SlipCosts {
    fn from_model(model: &StutterModel) -> Self {
        let ln_same_length_share = model.same_length_share().ln();
        Self {
            open_expansion: (model.whole_repeat_longer_share()
                * model.whole_repeat_one_step_share())
            .ln()
                - ln_same_length_share,
            open_contraction: (model.whole_repeat_shorter_share()
                * model.whole_repeat_one_step_share())
            .ln()
                - ln_same_length_share,
            extend: (1.0 - model.whole_repeat_one_step_share()).ln(),
        }
    }
}

/// Per-worker scratch: a ring of `period + 1` score rows (a whole-unit insertion reaches
/// `(i, j)` from `(i − period, j)`, so two rolling rows do not suffice) plus the full
/// backpointer matrix.
///
/// Grow-and-keep; buffers only, deciding nothing that changes a result.
#[derive(Debug, Default)]
pub struct AnchorRobustScratch {
    /// A ring of score rows, each `reference_len + 1` cells of `[f64; STATES]`.
    rows: Vec<Vec<[f64; STATES]>>,
    /// The winning predecessor state per cell and per state — flat, row-major.
    backpointers: Vec<[State; STATES]>,
}

impl AnchorRobustScratch {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Size the ring and the backpointer matrix. Grows only.
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

/// Pick the best of the candidates, keeping the **first on ties** — candidates are passed in
/// priority order, so the caller encodes the tie-break by their order.
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

/// Algorithm 5r: algorithm 5's firm anchor test over a boundary that a single miscall cannot
/// move. Generic over its [`Emission`], never behind `dyn`.
#[derive(Debug, Clone, Copy)]
pub struct SsrAnchorRobustAligner<E: Emission> {
    emission: E,
    costs: TransitionCosts,
    guard: JunctionGuard,
    anchor: AnchorRule,
}

impl<E: Emission> SsrAnchorRobustAligner<E> {
    /// The shipped configuration — [`JunctionGuard::default`] and [`AnchorRule::default`].
    #[must_use]
    pub fn new(emission: E) -> Self {
        Self::with_config(emission, JunctionGuard::default(), AnchorRule::default())
    }

    /// An explicit configuration, for experiments that need to report the settings they ran at.
    #[must_use]
    pub fn with_config(emission: E, guard: JunctionGuard, anchor: AnchorRule) -> Self {
        Self {
            emission,
            costs: TransitionCosts::new(),
            guard,
            anchor,
        }
    }

    /// Measure the read's repeat.
    ///
    /// `None` when there is no reference frame — a defined answer, not an error.
    #[must_use]
    pub fn delimit(
        &self,
        read: ReadBases<'_>,
        reference: &[u8],
        context: &RepeatContext<'_>,
        scratch: &mut AnchorRobustScratch,
    ) -> Option<super::ssr_best_path_flat_gap::TractReadout> {
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

        let slip = SlipCosts::from_model(context.stutter);
        let ring_len = period + 1;
        scratch.resize(read_len, reference_len, ring_len);
        let stride = reference_len + 1;
        let insertion_emission = self.emission.insert_ln();

        // A reference column is inside the tract when a gap touching it is a tract gap.
        let column_in_tract = |column: usize| left_flank_len < column && column <= tract_last;

        let guard = self.guard.resolve(
            tract_last.saturating_sub(left_flank_len),
            left_flank_len,
            right_flank_len,
        );

        // **The boundary guard.** The per-base gap-open cost by column: the ordinary cost away
        // from the junctions, and `UNREACHABLE` in the guarded columns either side of each one.
        // Flank columns outside the guard keep the ordinary open cost, which is what leaves a
        // genuine flank indel representable.
        let gap_open = |column: usize| {
            if column_in_tract(column) {
                let near_junction =
                    column <= left_flank_len + guard.tract || column + guard.tract > tract_last;
                if guard.tract > 0 && near_junction {
                    UNREACHABLE
                } else {
                    self.costs.ln_gap_open_tract()
                }
            } else {
                // Each test must name **which** flank it is about. A left-flank column is one at
                // or before the left junction; a right-flank column is one past the tract's last.
                // Testing only the distance would make every left-flank column "near" the right
                // junction as well, which bans the whole flank and throws away the flank-indel
                // route this algorithm exists to keep.
                let near_left = guard.left_flank > 0
                    && column <= left_flank_len
                    && column + guard.left_flank > left_flank_len;
                let near_right = guard.right_flank > 0
                    && column > tract_last
                    && column <= tract_last + guard.right_flank;
                if near_left || near_right {
                    UNREACHABLE
                } else {
                    self.costs.ln_gap_open()
                }
            }
        };

        // The **terminal** deletion route, exempt from the guard: a read that stops inside the
        // tract still has to reach the frame's end, and one that starts inside it still has to
        // be reached from the frame's start. Those are the `FromLeft` / `FromRight` cases, and
        // closing them would turn every partial read into a wrong answer.
        let terminal_gap_open = self.costs.ln_gap_open();

        // The `k`-th base of a whole unit inserted at tract column `column`: the motif,
        // phase-aligned so a unit inserted at a unit boundary begins at `motif[0]`.
        let motif_base = |column: usize, k: usize| motif[((column - left_flank_len) + k) % period];

        let AnchorRobustScratch { rows, backpointers } = scratch;

        // Row 0 — no read base consumed. Every deletion here is a *leading* one (the read starts
        // later in the frame), so it takes the terminal cost regardless of the guard.
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
                // Whole-unit deletion — only when the deleted unit lies wholly in the tract.
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

        // Rows 1..=read_len.
        for row_index in 1..=read_len {
            let read_base = bases[row_index - 1];
            let scores = self.emission.scores_for(read.quality_at(row_index - 1));
            let back_row = row_index * stride;

            let prev_slot = (row_index - 1) % ring_len;
            let slip_slot = row_index.checked_sub(period).map(|r| r % ring_len);
            let cur_slot = row_index % ring_len;
            // A deletion on the last read row is a *trailing* one — the read has run out and the
            // rest of the frame must still be crossed. That is the `FromLeft` case, not the
            // boundary-sliding gap the guard is aimed at.
            let trailing = row_index == read_len;

            // Column 0 — a read base before any reference base: the leading-insertion route,
            // which is how a read overhanging the frame's start is placed. Kept at the ordinary
            // open cost for the same reason the terminal deletions are.
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
                    // A completed slip closes to a match at no extra cost — the geometric has
                    // already priced the whole run.
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

                // Out-of-frame single-base insertion (a read base, no reference base).
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

                // Out-of-frame single-base deletion (a reference base, no read base). On the last
                // read row the terminal route applies, so a guarded open reopens.
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
                // motif, at a tract column. Each inserted base is scored at **its own** quality.
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
                // unit lies wholly in the tract. No emissions.
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

        // Final cell (m, n): best terminal state. Tie-break order M > D > I > SlipIns > SlipDel.
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
        };
        Some(trace_back(
            final_state,
            backpointers,
            geometry,
            Aligned {
                read: bases,
                reference,
            },
            self.anchor,
        ))
    }
}

/// The matrix shape the traceback needs — bundled so the walk takes one `Copy` argument instead
/// of six interchangeable `usize`s.
#[derive(Debug, Clone, Copy)]
struct MatrixGeometry {
    stride: usize,
    read_len: usize,
    reference_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
    period: usize,
}

/// The two sequences the traceback compares, bundled so the anchor test can ask whether a
/// `Match` step really matched.
#[derive(Debug, Clone, Copy)]
struct Aligned<'a> {
    read: &'a [u8],
    reference: &'a [u8],
}

impl Aligned<'_> {
    /// Whether read base `read_index` and reference base `reference_index` are the *same* base —
    /// the evidence unit the anchor test counts.
    ///
    /// Case-insensitive, because neither a read nor a reference is guaranteed upper-case. An `N`
    /// on either side is never evidence: it is the absence of a base call.
    #[inline]
    fn bases_agree(&self, read_index: usize, reference_index: usize) -> bool {
        let read = self.read[read_index].to_ascii_uppercase();
        let reference = self.reference[reference_index].to_ascii_uppercase();
        read == reference && read != b'N'
    }
}

/// The evidence a flank must show before it counts as anchoring: how far from the junction to
/// look, and how many of those bases must have *agreed*.
///
/// **Constructor state, not scratch** — both change results.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnchorRule {
    /// How many reference bases either side of the tract the test looks at — the window abutting
    /// the junction, where the bases that pin a boundary live.
    pub window: usize,
    /// How many bases in that window must agree. A quorum, not the whole window: a confident
    /// miscall at the junction must not demote a genuinely spanning read. Capped per locus by
    /// the flank actually present.
    pub min_matches: usize,
}

impl Default for AnchorRule {
    /// **A majority of a one-unit window**: seven bases, four of which must agree.
    ///
    /// The arithmetic is what sets the quorum, and it is the single biggest term in this
    /// algorithm's score. A *fabricated* crossing — a read that ran out inside a long tract and
    /// spent its repeat bases against the flank — agrees with the flank only by chance, 1 base
    /// in 4. Four or more agreements in seven such bases happens 7% of the time; algorithm 5's
    /// two-in-five happens 37% of the time. So five sixths of the fabricated completes that
    /// algorithm 5 waves through are caught here. A *genuine* spanning read at a 4% miscall rate
    /// expects 0.3 errors across the window and clears four essentially always, which is why
    /// the tightening costs the honest read nothing.
    ///
    /// The ceiling is set by how much flank a read has to show at all: production admits a
    /// complete read on 5 bp of flank, so a quorum above five would make those reads
    /// unmeasurable no matter how clean they are — and the bake-off confirms it, with the clean
    /// axis collapsing to 0.74 at a quorum of six. Four leaves those reads one miscall of slack.
    fn default() -> Self {
        Self {
            window: 7,
            min_matches: 4,
        }
    }
}

/// What the traceback found at each junction, and the rule that turns it into a classification.
///
/// Two questions per side, deliberately kept apart: did the path **cross** the junction
/// (positional, algorithm 4's test), and did the read **match the flank** beside it (evidential,
/// algorithm 5's test). Algorithm 5 required both on each side independently; the difference
/// here is entirely in [`Anchoring::decide`].
#[derive(Debug, Clone, Copy)]
struct Anchoring {
    left_crossed: bool,
    right_crossed: bool,
    left_matches: usize,
    right_matches: usize,
    left_needed: usize,
    right_needed: usize,
}

impl Anchoring {
    /// The classification, as a floor rather than as two independent tests.
    ///
    /// A crossing that the quorum rejects is dropped **only when something is left standing**.
    /// The three cases, in the order they are decided:
    ///
    /// 1. **The path crossed one junction only.** Whatever the quorum says, the answer is
    ///    already a lower bound, and the only demotion available is to nothing at all. A read
    ///    that reached a junction has located itself in the frame; a bound is a usable
    ///    observation and `Unanchored` is not. The crossing stands.
    /// 2. **The path crossed both, and the quorum accepts at least one.** The accepted sides
    ///    stand and the rejected one is dropped — a complete claim demoted to a bound, which is
    ///    exactly the over-claim algorithm 5 exists to catch.
    /// 3. **The path crossed both and the quorum accepts neither.** The evidence is too thin to
    ///    certify a length, but not a reason to discard the read: the better-evidenced side is
    ///    kept as a bound (the left on a tie, matching the M-first tie-break used throughout).
    ///
    /// The invariant across all three: this function never returns `(false, false)` for a read
    /// whose path crossed a junction. It can turn a measurement into a bound; it can never turn
    /// a bound into nothing.
    fn decide(self) -> (bool, bool) {
        let left_firm = self.left_crossed && self.left_matches >= self.left_needed;
        let right_firm = self.right_crossed && self.right_matches >= self.right_needed;

        match (self.left_crossed, self.right_crossed) {
            // Case 1 — at most a bound already; never demote further.
            (true, false) | (false, true) => (self.left_crossed, self.right_crossed),
            // Cases 2 and 3 — a complete claim, which the quorum may cut back to a bound.
            (true, true) => {
                if left_firm || right_firm {
                    (left_firm, right_firm)
                } else if self.left_matches >= self.right_matches {
                    (true, false)
                } else {
                    (false, true)
                }
            }
            (false, false) => (false, false),
        }
    }
}

/// Walk the traceback, reading the tract off the two flank junctions while counting, per flank,
/// the *agreeing* bases in the window abutting the junction.
fn trace_back(
    final_state: State,
    backpointers: &[[State; STATES]],
    geometry: MatrixGeometry,
    aligned: Aligned<'_>,
    anchor: AnchorRule,
) -> super::ssr_best_path_flat_gap::TractReadout {
    let MatrixGeometry {
        stride,
        read_len,
        reference_len,
        left_flank_len,
        right_flank_len,
        period,
    } = geometry;
    let left_junction = left_flank_len; // first tract reference base
    let right_junction = reference_len - right_flank_len; // first right-flank base
    let mut tract_start = 0usize;
    let mut tract_end = read_len;

    let left_window = left_junction.saturating_sub(anchor.window)..left_junction;
    let right_window = right_junction..(right_junction + anchor.window).min(reference_len);
    let mut left_matches = 0usize;
    let mut right_matches = 0usize;

    let mut i = read_len;
    let mut j = reference_len;
    let mut state = final_state;
    while i != 0 || j != 0 {
        let pred = backpointers[i * stride + j][state.index()];
        match state {
            State::Match => {
                let consumed = j - 1;
                if consumed == left_junction {
                    tract_start = i - 1;
                }
                if right_flank_len > 0 && consumed == right_junction {
                    tract_end = i - 1;
                }
                if (left_window.contains(&consumed) || right_window.contains(&consumed))
                    && aligned.bases_agree(i - 1, consumed)
                {
                    if consumed < left_junction {
                        left_matches += 1;
                    } else {
                        right_matches += 1;
                    }
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
                // A contraction at the very start of the tract deletes the first tract base,
                // which *is* the left junction, so the crossing must be recorded exactly as an
                // ordinary deletion records it.
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

    // The quorum is capped by the flank the locus actually has: a 1 bp flank can never supply
    // two agreeing bases, and demanding it would make every contig-edge locus unmeasurable
    // rather than merely thinly evidenced.
    let anchoring = Anchoring {
        left_crossed: left_flank_len > 0 && tract_start != 0,
        right_crossed: right_flank_len > 0 && tract_end != read_len,
        left_matches,
        right_matches,
        left_needed: anchor.min_matches.min(left_flank_len),
        right_needed: anchor.min_matches.min(right_flank_len),
    };
    let (left_anchored, right_anchored) = anchoring.decide();

    // A dropped junction leaves a *stale offset* behind, and the offset is load-bearing: a
    // `FromLeft` span must run to the end of the read, because that is what makes it a lower
    // bound on the allele. Keeping a crossing the quorum just rejected would report a bound
    // shorter than the repeat the read actually showed.
    if !left_anchored {
        tract_start = 0;
    }
    if !right_anchored {
        tract_end = read_len;
    }

    super::ssr_best_path_flat_gap::TractReadout {
        tract_start: tract_start as u64,
        tract_end: tract_end as u64,
        left_anchored,
        right_anchored,
    }
}

impl<E: Emission> BestPathAligner for SsrAnchorRobustAligner<E> {
    type Scratch = AnchorRobustScratch;
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
    use crate::ng::alignment::stutter::StutterRates;
    use crate::ng::alignment::{RepeatGeometry, StutterModel};
    use crate::ng::types::{Bp, Motif};

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

    /// Contraction-biased parameters — HipSTR's fitted values are.
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

    fn measure(
        read: &[u8],
        reference: &[u8],
        geometry: &RepeatGeometry,
        model: &StutterModel,
    ) -> RepeatSpan {
        let aligner = SsrAnchorRobustAligner::new(PerQualityEmission::new());
        let context = RepeatContext {
            geometry,
            stutter: model,
        };
        let quality = vec![35u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = AnchorRobustScratch::new();
        aligner.align(bases, reference, context, &mut scratch)
    }

    /// **A contraction must be cheaper to open than an expansion**, on a model that says so.
    /// The asymmetry lives at the cost level, and this is the only thing in this file that
    /// reads it: without it, `open_expansion` and `open_contraction` could be exchanged and
    /// every other test here would still pass.
    #[test]
    fn a_contraction_is_cheaper_to_open_than_an_expansion() {
        let slip = SlipCosts::from_model(&contraction_biased());
        assert!(
            slip.open_contraction > slip.open_expansion,
            "a contraction-biased model must make contraction the cheaper slip to open"
        );
    }

    /// The slip costs reconstruct the stutter model's own probabilities exactly — the affine
    /// open plus its extends must sum back to `ln(P(n)) − ln(same_length_share)`, or this
    /// aligner is pricing a different distribution than the one it was handed.
    #[test]
    fn the_slip_costs_reconstruct_the_stutter_probability() {
        let model = contraction_biased();
        let slip = SlipCosts::from_model(&model);
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

    /// **Mandatory property 1: a clean read of the reference measures the reference length.**
    #[test]
    fn a_clean_read_measures_the_reference_tract() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let span = measure(&reference, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
    }

    /// A genuine whole-unit expansion is measured at its own length — the slip route must stay
    /// usable through the guard, which only ever closes the *out-of-frame* route.
    #[test]
    fn a_whole_unit_expansion_is_measured() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAGCAGCAGCAGTTGGTTGGAT"; // +2 units
        let span = measure(read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(18));
    }

    /// **Mandatory property 2: out-of-frame changes keep a route.** An interrupted repeat must
    /// still be measured verbatim, not forced onto a tidy unit count — the clamp on the tract
    /// guard exists precisely so this stays true.
    #[test]
    fn an_out_of_frame_change_still_has_a_route() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCTGCAGTTGGTTGGAT";
        let span = measure(read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
        let observed = span.observed_span().expect("a measured span");
        assert_eq!(
            &read[observed.start as usize..observed.end as usize],
            b"CAGCAGCTGCAG"
        );
    }

    /// **A one-base out-of-frame insertion in the tract interior survives the guard.** This is
    /// the route the tract guard closes near the junctions, and it must remain open away from
    /// them — otherwise an impure repeat is rounded onto a whole unit count.
    #[test]
    fn an_interior_one_base_insertion_survives_the_guard() {
        let motif = b"AAG";
        let left = b"GATCTTGCAAGCTGGAATCCGTTAC";
        let right = b"CAGTTCACGATCCTAAGGCTTGACG";
        let ref_tract: Vec<u8> = motif.iter().cycle().take(24).copied().collect();
        let (reference, geometry) = frame(left, &ref_tract, right, motif);

        let mut sample = ref_tract.clone();
        sample.insert(12, b'T'); // mid-tract, well clear of both guards
        let mut read = left.to_vec();
        read.extend_from_slice(&sample);
        read.extend_from_slice(right);

        let span = measure(&read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(25));
    }

    /// **The flank-indel route stays open — the reason this is not algorithm 4n.** A real 1 bp
    /// indel in the body of a flank must still be spellable, so the tract measurement is
    /// unaffected by it. Both sides, insertion and deletion.
    #[test]
    fn a_real_flank_indel_does_not_derail_the_measurement() {
        let motif = b"AAG";
        let left = b"GATCTTGCAAGCTGGAATCCGTTAC";
        let right = b"CAGTTCACGATCCTAAGGCTTGACG";
        let ref_tract: Vec<u8> = motif.iter().cycle().take(24).copied().collect();
        let (reference, geometry) = frame(left, &ref_tract, right, motif);

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
                let (read_left, read_right) = if right_side {
                    (left.to_vec(), edit(right))
                } else {
                    (edit(left), right.to_vec())
                };
                let mut read = read_left;
                read.extend_from_slice(&ref_tract);
                read.extend_from_slice(&read_right);

                let span = measure(&read, &reference, &geometry, &contraction_biased());
                assert_eq!(
                    span.measured_length(),
                    Some(24),
                    "a flank indel derailed the measurement (right_side {right_side}, insert {insert})"
                );
            }
        }
    }

    /// **The reason the anchor test exists.** A long allele that outruns the read must come back
    /// as a lower bound, not a fabricated complete length.
    #[test]
    fn a_long_allele_that_outruns_the_read_is_a_lower_bound_not_a_measurement() {
        let motif = b"AAG";
        let left = b"GATCTTGCAAGCTGGAATCCGTTAC";
        let right = b"CAGTTCACGATCCTAAGGCTTGACG";
        let ref_tract: Vec<u8> = motif.iter().cycle().take(24).copied().collect();
        let (reference, geometry) = frame(left, &ref_tract, right, motif);

        let mut read = left.to_vec();
        read.extend(motif.iter().cycle().take(55).copied());

        let span = measure(&read, &reference, &geometry, &contraction_biased());
        assert!(
            matches!(span, RepeatSpan::FromLeft(_)),
            "a read that never reached the right flank must not report a length, got {span:?}"
        );
        assert_eq!(span.measured_length(), None);
        let observed = span
            .observed_span()
            .expect("a lower bound still has a span");
        assert_eq!(observed.end, read.len() as u64);
    }

    /// The mirror case: a read that begins inside a long tract and runs into the right flank.
    #[test]
    fn a_read_starting_inside_a_long_tract_anchors_only_on_the_right() {
        let motif = b"AAG";
        let left = b"GATCTTGCAAGCTGGAATCCGTTAC";
        let right = b"CAGTTCACGATCCTAAGGCTTGACG";
        let ref_tract: Vec<u8> = motif.iter().cycle().take(24).copied().collect();
        let (reference, geometry) = frame(left, &ref_tract, right, motif);

        let mut read: Vec<u8> = motif.iter().cycle().take(55).copied().collect();
        read.extend_from_slice(right);

        let span = measure(&read, &reference, &geometry, &contraction_biased());
        assert!(
            matches!(span, RepeatSpan::FromRight(_)),
            "a read that never reached the left flank must not report a length, got {span:?}"
        );
        let observed = span
            .observed_span()
            .expect("a lower bound still has a span");
        assert_eq!(observed.start, 0);
    }

    /// **The floor: a bound is never demoted to nothing.** This is the one classification rule
    /// that differs from algorithm 5, and it is checked at the decision itself because the
    /// combination it guards against (a crossing with no agreeing flank bases) is reachable
    /// under noise but awkward to provoke deterministically through the whole aligner.
    #[test]
    fn a_crossed_junction_is_never_demoted_to_nothing() {
        for left_crossed in [false, true] {
            for right_crossed in [false, true] {
                for left_matches in 0..=3 {
                    for right_matches in 0..=3 {
                        let anchoring = Anchoring {
                            left_crossed,
                            right_crossed,
                            left_matches,
                            right_matches,
                            left_needed: 2,
                            right_needed: 2,
                        };
                        let (l, r) = anchoring.decide();
                        if left_crossed || right_crossed {
                            assert!(
                                l || r,
                                "a crossed junction was demoted to nothing: {anchoring:?}"
                            );
                        }
                        // And nothing is ever *promoted*: an uncrossed junction stays uncrossed.
                        assert!(!l || left_crossed, "the left junction was promoted");
                        assert!(!r || right_crossed, "the right junction was promoted");
                    }
                }
            }
        }
    }

    /// A single crossing keeps its side no matter how thin the flank evidence is — the exact
    /// case that made algorithm 5 lose real partial reads.
    #[test]
    fn a_lone_crossing_survives_a_failed_quorum() {
        let anchoring = Anchoring {
            left_crossed: true,
            right_crossed: false,
            left_matches: 0,
            right_matches: 0,
            left_needed: 2,
            right_needed: 2,
        };
        assert_eq!(anchoring.decide(), (true, false));
    }

    /// A complete claim with evidence on one side only is cut back to that side — the demotion
    /// the anchor test exists for still fires.
    #[test]
    fn a_one_sided_quorum_demotes_a_complete_claim_to_a_bound() {
        let anchoring = Anchoring {
            left_crossed: true,
            right_crossed: true,
            left_matches: 5,
            right_matches: 0,
            left_needed: 2,
            right_needed: 2,
        };
        assert_eq!(anchoring.decide(), (true, false));
    }

    /// **The demotion must not fire on real reads.** A genuinely spanning read with a
    /// sequencing error in each junction window still measures its allele.
    #[test]
    fn a_miscall_beside_each_junction_does_not_demote_a_spanning_read() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let mut read = reference.clone();
        read[7] = b'A'; // last base of the left flank, miscalled
        read[20] = b'A'; // first base of the right flank, miscalled
        let span = measure(&read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
    }

    /// **The noise fix, at the mechanism.** A miscall immediately inside each junction is
    /// exactly the configuration that buys a one-base slide when a gap may open beside it. With
    /// the guard the boundary holds and the length is right; the assertion is the length,
    /// because the length is what the caller consumes.
    #[test]
    fn a_miscall_inside_each_junction_does_not_move_the_boundary() {
        let motif = b"AAG";
        let left = b"GATCTTGCAAGCTGGAATCCGTTAC";
        let right = b"CAGTTCACGATCCTAAGGCTTGACG";
        let ref_tract: Vec<u8> = motif.iter().cycle().take(24).copied().collect();
        let (reference, geometry) = frame(left, &ref_tract, right, motif);

        let mut read = reference.to_vec();
        let tract_start = left.len();
        let tract_end = left.len() + ref_tract.len();
        read[tract_start] = b'C'; // first tract base, miscalled
        read[tract_end - 1] = b'C'; // last tract base, miscalled

        let span = measure(&read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(24));
    }

    /// A locus whose flanks are shorter than the quorum is judged against the flank it has.
    #[test]
    fn a_flank_shorter_than_the_quorum_can_still_anchor() {
        let (reference, geometry) = frame(b"T", b"CAGCAGCAGCAG", b"T", b"CAG");
        let span = measure(&reference, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
    }

    /// **The clamp is load-bearing.** On a short tract or a stub flank an unclamped guard would
    /// close the whole out-of-frame route. Configuring an absurd guard must therefore change
    /// nothing that the clamp can prevent: the resolved widths always leave interior behind.
    #[test]
    fn the_guard_is_clamped_so_the_tract_and_flanks_keep_a_route() {
        for tract_len in [3usize, 6, 9, 12, 24, 60] {
            for flank_len in [0usize, 1, 3, 5, 25] {
                let resolved = JunctionGuard {
                    tract_bases: 100,
                    flank_bases: 100,
                }
                .resolve(tract_len, flank_len, flank_len);
                assert!(
                    2 * resolved.tract < tract_len || resolved.tract == 0,
                    "the tract guard closed the whole interior at tract_len {tract_len}"
                );
                assert!(
                    2 * resolved.left_flank <= flank_len,
                    "the flank guard closed the whole left flank at flank_len {flank_len}"
                );
                assert!(2 * resolved.right_flank <= flank_len);
            }
        }
    }

    /// The shipped flank guard must stay clear of where a real flank indel sits. Production
    /// requires 5 bp of flank on a complete read, and the bake-off plants its indel in the
    /// middle of a 25 bp flank; a guard that reached that far would trade the flank-indel axis
    /// for the noise axis, which is the trade this algorithm exists to avoid.
    #[test]
    fn the_shipped_flank_guard_leaves_the_flank_body_open() {
        let resolved = JunctionGuard::default().resolve(24, 25, 25);
        assert!(
            resolved.left_flank < 25 / 2,
            "the flank guard reaches the flank body, where real indels live"
        );
        assert_eq!(resolved.left_flank, resolved.right_flank);
    }

    /// The shipped quorum must fit inside the window it is counted in — a quorum larger than its
    /// window can never be met, so every flank would read as unanchored — and a quorum of zero
    /// is no test at all.
    #[test]
    fn the_anchor_quorum_fits_inside_the_anchor_window() {
        let rule = AnchorRule::default();
        assert!(
            rule.min_matches <= rule.window,
            "a quorum larger than the window it is counted in can never be met"
        );
        assert!(rule.min_matches > 0, "a zero quorum is no test at all");
    }

    /// An empty reference is a defined answer, not an error.
    #[test]
    fn an_empty_reference_is_unanchored() {
        let geometry = RepeatGeometry {
            left_flank_len: Bp(0),
            right_flank_len: Bp(0),
            motif: Motif::new(b"CAG").expect("valid"),
        };
        assert_eq!(
            measure(b"ACGT", b"", &geometry, &contraction_biased()),
            RepeatSpan::Unanchored
        );
    }

    /// Scratch reuse across a size and period drop must not leak.
    #[test]
    fn scratch_reuse_does_not_leak_across_periods() {
        let aligner = SsrAnchorRobustAligner::new(PerQualityEmission::new());
        let model = contraction_biased();

        let loci: &[(&[u8], &[u8])] = &[
            (b"GTTGTG", b"GTTGTGGTTGTGGTTGTGGTTGTGGTTGTG"),
            (b"A", b"AAAA"),
            (b"CAG", b"CAGCAGCAGCAG"),
            (b"CA", b"CACACACACACA"),
        ];

        let mut reused_scratch = AnchorRobustScratch::new();
        for &(motif, tract) in loci {
            let (reference, geometry) = frame(b"ACGTACGT", tract, b"TTGGTTGGAT", motif);
            let quality = vec![35u8; reference.len()];
            let context = RepeatContext {
                geometry: &geometry,
                stutter: &model,
            };
            let bases = ReadBases::try_new(&reference, &quality).expect("matched lengths");
            let reused = aligner.align(bases, &reference, context, &mut reused_scratch);
            let bases = ReadBases::try_new(&reference, &quality).expect("matched lengths");
            let fresh = aligner.align(bases, &reference, context, &mut AnchorRobustScratch::new());
            assert_eq!(
                reused,
                fresh,
                "scratch reuse leaked across periods at motif {:?}",
                std::str::from_utf8(motif).unwrap()
            );
        }
    }

    /// Multiple periods work — period 1 (homopolymer) and period 2, where the slip stride is
    /// shortest and the guard is widest relative to the tract.
    #[test]
    fn homopolymer_and_dinucleotide_loci_are_measured() {
        let (reference, geometry) = frame(b"ACGTACGT", b"AAAAAAAAAAAA", b"TTGGTTGGAT", b"A");
        assert_eq!(
            measure(&reference, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(12)
        );
        let read = b"ACGTACGTAAAAAAAAAAAAATTGGTTGGAT"; // +1 base
        assert_eq!(
            measure(read, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(13)
        );

        let (reference, geometry) = frame(b"ACGTACGT", b"CACACACACACA", b"TTGGTTGGAT", b"CA");
        assert_eq!(
            measure(&reference, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(12)
        );
    }

    /// **The load-bearing discriminants.** `State` writes score arrays by position and reads them
    /// by `index()`, so reordering the variants would silently permute every cell.
    #[test]
    fn the_state_discriminants_are_the_array_indices() {
        assert_eq!(State::Match.index(), 0);
        assert_eq!(State::Insertion.index(), 1);
        assert_eq!(State::Deletion.index(), 2);
        assert_eq!(State::SlipInsertion.index(), 3);
        assert_eq!(State::SlipDeletion.index(), 4);
    }
}
