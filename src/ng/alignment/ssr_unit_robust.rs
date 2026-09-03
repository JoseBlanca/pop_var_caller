//! Algorithm 4r — **unit-robust**: algorithm 4's recurrence and scoring, with two changes that
//! attack the two ways a delimiter gets a tract boundary wrong, and **without banning flank
//! gaps**.
//!
//! The scoring model is algorithm 4's, unchanged: a whole-unit slip priced from the shared
//! [`StutterModel`] plus a flat per-base gap for everything out of frame (see
//! [`super::ssr_best_path_unit_slip`] for the derivation). What changes is *where a per-base gap
//! may open* and *when a crossed junction counts as an anchor*.
//!
//! # The two failures, and the two answers
//!
//! **1. The fabricated complete.** Algorithm 4 calls a flank anchored from the *shape* of the
//! traceback alone: the left flank held if the tract did not start at read offset 0, the right
//! flank held if it did not end at the read's last base. That asks whether the path *crossed* the
//! junction, not whether the read had any *evidence* for the flank it crossed into. On a long
//! allele that outruns the read there is no right flank in the read at all, yet the frame still
//! has 25 flank bases to account for: spending the read's last repeat bases as mismatched
//! "matches" against them is cheaper than a 25-base deletion run. The path crosses, algorithm 4
//! says `Between`, and reports a **length for an allele the read never finished reading**. The
//! answer is [`the firm-anchor test`](trace_back): a flank anchors only when the read *matched
//! its bases* in the window abutting the junction ([`AnchorRule`]).
//!
//! **2. The one-mismatch boundary slide.** At Q40 a mismatch scores `ln(1e-4/3) ≈ -10.31` and the
//! flank gap-open scores `ln(2.9e-5) ≈ -10.45`. **One bad base costs about what opening a gap
//! costs**, so an insertion plus a deletion — which together slide the read one base against the
//! reference — is worth taking whenever the slide repairs three miscalls. The slide moves the
//! junction and the measurement is wrong by a base or, at period 1, by a whole unit. The answer is
//! the [`JunctionGuard`]: within a few columns *either side* of each junction, the out-of-frame
//! per-base route may not **open**. Nothing else changes — it may still extend, it may still open
//! anywhere else, and whole-unit slips are never touched, so genuine length variation is
//! untouched.
//!
//! # Why the guard is a window and not a flank ban
//!
//! Algorithm 4n ([`super::ssr_noise_robust`]) closes the same door by banning per-base gaps in the
//! flank outright. That works on noise and **destroys real flank indels**: a read whose flank
//! genuinely carries a 1 bp indel has no route to express it, so the gap reappears inside the
//! tract and the tract is mismeasured — 4n scores 0.000 on the flank-indel axis of the synthetic
//! bake-off. The failure is not the guarding, it is the *width*: what pins a boundary is the
//! handful of bases beside it, and an indel 12 bases into a flank is not a boundary event. So the
//! guard here is a **window around each junction**, on both the flank side and the tract side, and
//! the rest of the flank keeps the ordinary gap-open. Real flank indels stay expressible; the
//! one-base slide at the junction does not.
//!
//! # Demotion is one-way, and it never removes the last anchor
//!
//! Both changes can only ever make a span **less** of a claim, never more, and this algorithm goes
//! one step further than algorithm 5 ([`super::ssr_anchor_firm`]): a read that algorithm 4 would
//! report as one-sided is never demoted at all, and a `Between` is never demoted on *both* sides.
//! The rule is that a **partial must not become [`RepeatSpan::Unanchored`]**. A lower bound is a
//! usable fact — the allele is at least this long — while `Unanchored` is the read thrown away, so
//! turning a partial into nothing is not a conservative correction, it is a lost read. The one
//! over-claim worth undoing is the fabricated complete, and that is exactly `Between` → one-sided.
//! See [`resolve_anchors`].

use super::emission::Emission;
use super::ssr_best_path_flat_gap::{TractReadout, TransitionCosts};
use super::stutter::StutterModel;
use super::{BestPathAligner, ReadBases, RepeatContext, RepeatSpan};
use crate::ng::types::MAX_MOTIF_LEN;

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

    /// The state a 3-bit field holds. Total over the field's domain, because a packed
    /// backpointer is read back with no other check — an out-of-range code would mean the
    /// packing is broken, and answering `Match` there would hide it behind a plausible
    /// traceback, so it panics instead.
    #[inline]
    const fn from_code(code: u16) -> Self {
        match code {
            0 => State::Match,
            1 => State::Insertion,
            2 => State::Deletion,
            3 => State::SlipInsertion,
            4 => State::SlipDeletion,
            _ => panic!("backpointer field out of range — the packing is broken"),
        }
    }
}

/// Bits per state in a packed backpointer cell. Five states need three.
const STATE_BITS: u32 = 3;

/// The five predecessors of one cell, three bits each, in one `u16`.
///
/// The DP writes one of these per cell and the traceback reads at most `read_len +
/// reference_len` of them, so the array is written orders of magnitude more often than it is
/// read: at the measured locus dimensions it is ~10⁴ cells per read and ~10⁷ cells per deep
/// locus. As `[State; 5]` that was five separate `strb` per cell — the record is 5 bytes, so
/// it is never naturally aligned and the stores never merge — and a `×5` address multiply.
/// Packed it is one `strh` at a power-of-two stride.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct Backpointers(u16);

// The packing's one invariant: five 3-bit fields must fit the `u16`, and every state code must
// fit its field. A sixth state would silently truncate, so it fails the build instead.
const _: () = assert!(STATES as u32 * STATE_BITS <= u16::BITS);
const _: () = assert!((STATES as u16) <= 1 << STATE_BITS);

impl Backpointers {
    /// Pack the five predecessors, in state order.
    #[inline]
    fn new(predecessors: [State; STATES]) -> Self {
        let mut packed = 0u16;
        for (index, state) in predecessors.iter().enumerate() {
            packed |= (*state as u16) << (STATE_BITS * index as u32);
        }
        Self(packed)
    }

    /// The predecessor recorded for entering the cell in `state`.
    #[inline]
    fn get(self, state: State) -> State {
        let shift = STATE_BITS * state.index() as u32;
        State::from_code((self.0 >> shift) & ((1 << STATE_BITS) - 1))
    }
}

/// An unreachable score. A real value: some cells cannot be entered in some states, and the
/// junction guard also uses it to *close* a route (a guarded gap-open is an unreachable
/// transition, not merely an expensive one).
const UNREACHABLE: f64 = f64::NEG_INFINITY;

/// How wide the no-open window around each junction is, in reference columns, on each side.
///
/// **Constructor state, not scratch**, so a bake-off can name the setting it compared at.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct JunctionGuard {
    /// Columns *inside the tract* at each end that may not open an out-of-frame gap.
    ///
    /// Clamped per locus so the tract keeps at least one motif unit of openable columns in the
    /// middle ([`Self::tract_width`]) — without that, a three-unit tract (the common case in a
    /// real genome) would have its whole interior closed and an interrupted repeat there would be
    /// forced onto a tidy unit count, which is the failure the out-of-frame route exists to avoid.
    pub tract_columns: usize,
    /// Columns *inside a flank*, abutting the junction, that may not open an out-of-frame gap.
    ///
    /// This is the knob that separates this algorithm from a flank ban. Clamped to half the flank
    /// ([`Self::flank_width`]) so a genuine indel in the body of the flank always keeps a route:
    /// the guard is for the bases that pin the boundary, not for the flank as a whole.
    pub flank_columns: usize,
}

impl JunctionGuard {
    /// The tract-side width actually usable at a tract of `tract_len` bases and motif `period` —
    /// the configured width cut back so at least one motif unit of interior columns can still
    /// open a gap.
    #[inline]
    #[must_use]
    fn tract_width(self, tract_len: usize, period: usize) -> usize {
        self.tract_columns.min(tract_len.saturating_sub(period) / 2)
    }

    /// The flank-side width actually usable at a flank of `flank_len` bases — never more than
    /// half of it, so the flank always keeps a body in which a real indel can open.
    #[inline]
    #[must_use]
    fn flank_width(self, flank_len: usize) -> usize {
        self.flank_columns.min(flank_len / 2)
    }
}

impl Default for JunctionGuard {
    /// Eight columns of tract and eight of flank on each side of each junction.
    ///
    /// The lower bound is principled: a sub-unit gap opening within one motif unit of the junction
    /// is a boundary convention rather than an interior event, and the widest [`Motif`] period is
    /// six, so six is the narrowest width that can be argued for. The upper bound is empirical —
    /// on the synthetic bake-off (`examples/ng_ssr_synthetic_bakeoff.rs`) the composite is **flat
    /// from six to twelve** (0.9993–0.9994) and only the clamps stop it drifting further. Eight
    /// sits in the middle of that plateau, which is the point least sensitive to being slightly
    /// wrong in either direction; taking the narrow edge would leave nothing in hand if a real
    /// locus behaves a little worse than the synthetic one.
    ///
    /// [`Motif`]: crate::ng::types::Motif
    fn default() -> Self {
        Self {
            tract_columns: 8,
            flank_columns: 8,
        }
    }
}

/// When a crossed junction counts as an anchor — the evidence test that turns a fabricated
/// complete back into the lower bound the read actually supports.
///
/// **Constructor state, not scratch**, for the same reason as [`JunctionGuard`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AnchorRule {
    /// How many reference bases either side of the tract the test looks at — the window abutting
    /// the junction, where the bases that pin a boundary live.
    ///
    /// Five by default, production's own admission floor for a complete read, so the test never
    /// demands evidence the upstream gate did not require.
    pub window: usize,
    /// How many bases in that window must **agree** for the flank to count as anchoring.
    ///
    /// A quorum, not the whole window: a single confident miscall at the junction must not demote
    /// a genuinely spanning read. Capped per locus by the flank that actually exists, so a repeat
    /// at a contig edge is judged against the flank it has rather than made unanchorable.
    ///
    /// Three of five is the balance point, and the two axes it balances pull in opposite
    /// directions: a quorum of two lets a noisy long allele buy a spurious anchor (the bake-off's
    /// `p_noise` falls to 0.988), while four starts demoting genuinely spanning noisy reads
    /// (`noise` falls from 0.999 to 0.990). Three loses neither.
    pub min_matches: usize,
    /// How many read bases must be **aligned against** a flank — anywhere in it, agreeing or not
    /// — before that flank counts as anchoring.
    ///
    /// The counterpart to [`Self::min_matches`], counted in columns rather than agreements: a
    /// heavily miscalled but genuinely present flank still anchors, while a flank represented by
    /// one lone stray base does not.
    ///
    /// Five, which is **exactly** production's admission floor for a complete read, and the value
    /// is a ceiling rather than a preference: six breaks the bake-off's clean axis outright
    /// (1.000 → 0.743), because a read showing only the five flank bases production admits can
    /// never supply six. So this rule is pinned to the upstream gate — it rejects the degenerate
    /// one- and two-base anchors and nothing else, and raising it would be relitigating where a
    /// read is admitted, from the wrong end of the pipeline.
    pub min_support: usize,
}

impl Default for AnchorRule {
    fn default() -> Self {
        Self {
            window: 5,
            min_matches: 3,
            min_support: 5,
        }
    }
}

/// The two knobs together — what separates algorithm 4r from algorithm 4.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct UnitRobustConfig {
    /// Where an out-of-frame gap may not open.
    pub guard: JunctionGuard,
    /// When a crossed junction counts as an anchor.
    pub anchor: AnchorRule,
}

/// The whole-unit slip transition costs, in log space, derived **from the shared
/// [`StutterModel`]** — algorithm 4's derivation, unchanged.
#[derive(Debug, Clone, Copy)]
struct SlipCosts {
    /// `ln(whole_repeat_longer_share · whole_repeat_one_step_share) − ln(same_length_share)`.
    open_expansion: f64,
    /// `ln(whole_repeat_shorter_share · whole_repeat_one_step_share) − ln(same_length_share)`.
    open_contraction: f64,
    /// `ln(1 − whole_repeat_one_step_share)`, per unit after the first.
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

/// What a reference column is worth, independent of the read — resolved once per call and read
/// once per cell (see the `plan` block in [`SsrUnitRobustAligner::delimit`]).
#[derive(Debug, Clone, Copy)]
struct ColumnPlan {
    /// The per-base gap-open cost at this column: tract-aware, and `UNREACHABLE` where the
    /// junction guard closes it.
    gap_open: f64,
    /// The same, reopened for a terminal route (`terminal_gap_open.max(gap_open)`) — a read that
    /// runs off an end must still cross the rest of the frame.
    gap_open_terminal: f64,
    /// Whether a gap touching this column is a tract gap.
    in_tract: bool,
    /// Whether a whole unit of reference bases ending at this column lies wholly in the tract —
    /// the precondition of the whole-unit deletion route.
    unit_deletable: bool,
}

/// Per-worker scratch: a ring of `period + 1` score rows, the full backpointer matrix (exactly as
/// algorithm 4 needs — a whole-unit insertion reaches `(i, j)` from `(i − period, j)`), and the
/// per-column plan.
///
/// Grow-and-keep; buffers only, deciding nothing that changes a result.
#[derive(Debug, Default)]
pub struct UnitRobustScratch {
    rows: Vec<Vec<[f64; STATES]>>,
    backpointers: Vec<Backpointers>,
    column_plan: Vec<ColumnPlan>,
}

impl UnitRobustScratch {
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
            self.backpointers.resize(cells, Backpointers::default());
        }
    }
}

/// Pick the best candidate, keeping the **first on ties** — the spec's determinism rule, encoded
/// by the caller's argument order (match, then deletion, then insertion, then the slips).
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

/// Algorithm 4r: the unit-robust repeat delimiter. Generic over its [`Emission`], never behind
/// `dyn` (arch §4).
#[derive(Debug, Clone, Copy)]
pub struct SsrUnitRobustAligner<E: Emission> {
    emission: E,
    costs: TransitionCosts,
    config: UnitRobustConfig,
}

impl<E: Emission> SsrUnitRobustAligner<E> {
    /// Build the delimiter at its shipped configuration.
    #[must_use]
    pub fn new(emission: E) -> Self {
        Self::with_config(emission, UnitRobustConfig::default())
    }

    /// Build the delimiter at an explicit configuration — the form a bake-off uses, so the
    /// setting under comparison is named at the call site.
    #[must_use]
    pub fn with_config(emission: E, config: UnitRobustConfig) -> Self {
        Self {
            emission,
            costs: TransitionCosts::new(),
            config,
        }
    }

    /// The configuration this delimiter is running at.
    #[must_use]
    pub fn config(&self) -> UnitRobustConfig {
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
        scratch: &mut UnitRobustScratch,
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

        let slip = SlipCosts::from_model(context.stutter);
        let ring_len = period + 1;
        scratch.resize(read_len, reference_len, ring_len);
        let stride = reference_len + 1;
        let insertion_emission = self.emission.insert_ln();

        // A reference column is inside the tract when a gap touching it is a tract gap.
        let column_in_tract = |column: usize| left_flank_len < column && column <= tract_last;

        // The junction guard, resolved to this locus: how many tract columns at each end, and how
        // many columns of each flank abutting its junction, may not *open* an out-of-frame gap.
        let guard = self.config.guard;
        let tract_guard = guard.tract_width(tract_last - left_flank_len, period);
        let left_guard = guard.flank_width(left_flank_len);
        let right_guard = guard.flank_width(right_flank_len);

        // Guarded columns, in the aligner's 1-based column space. The left junction sits between
        // columns `left_flank_len` and `left_flank_len + 1`; the right junction between
        // `tract_last` and `tract_last + 1`. Each guard reaches `width` columns either way.
        let guarded = |column: usize| {
            let near_left = column + left_guard > left_flank_len
                && column <= left_flank_len + tract_guard
                && (left_guard > 0 || tract_guard > 0);
            let near_right = column + tract_guard > tract_last
                && column <= tract_last + right_guard
                && (tract_guard > 0 || right_guard > 0);
            near_left || near_right
        };

        // The per-base gap-open cost by column: algorithm 4's tract-aware cost, closed at the
        // guarded columns. **Nothing else about the flank changes** — an indel in the body of a
        // flank keeps its ordinary open, which is what algorithm 4n gave up.
        let gap_open = |column: usize| {
            if guarded(column) {
                UNREACHABLE
            } else if column_in_tract(column) {
                self.costs.ln_gap_open_tract()
            } else {
                self.costs.ln_gap_open()
            }
        };

        // The **terminal** routes, exempt from the guard: a read that stops inside the tract still
        // has to reach the frame's end, and one that starts inside it still has to be reached from
        // the frame's start. Those crossings are what `FromLeft`/`FromRight` report, so closing
        // them would turn every partial read into a wrong answer. Recognised structurally — a
        // deletion on read row 0 or on the last read row — never by a threshold.
        let terminal_gap_open = self.costs.ln_gap_open();

        // Resolve the column axis **once per read** instead of once per cell.
        //
        // `gap_open`, `column_in_tract` and the whole-unit-deletion window are functions of the
        // column and the locus geometry alone — no read input — yet the cell body evaluated all
        // three, so each was re-derived `read_len` times per column: a `ccmp`/`cset`/`fcsel` chain
        // with four locus constants reloaded from the stack, plus an unrolled six-step compare
        // chain for the slip window. Here they cost one pass over the axis and the cell reads one
        // [`ColumnPlan`].
        let plan = &mut scratch.column_plan;
        plan.clear();
        plan.reserve(reference_len + 1);
        for column in 0..=reference_len {
            let open = gap_open(column);
            plan.push(ColumnPlan {
                gap_open: open,
                // The terminal routes reopen a guarded column, so both forms are wanted per cell.
                gap_open_terminal: terminal_gap_open.max(open),
                in_tract: column_in_tract(column),
                // A whole unit of reference bases can only be deleted when every column it covers
                // is a tract column — the range test the cell used to unroll per cell.
                unit_deletable: column >= period
                    && (column - period + 1..=column).all(column_in_tract),
            });
        }

        // The whole-unit slip emission, per read row, indexed by the motif **phase** the tract
        // column sits at — `(column - left_flank_len) % period`, so a unit inserted at a unit
        // boundary is scored against `motif[0]`.
        //
        // **Why a per-row table rather than a per-cell sum.** The sum runs over the `period` read
        // bases ending at this row, each scored against a motif base, so it depends on the row and
        // on the column's phase alone — and a phase takes only `period` values. Computed in the
        // cell it was re-derived at every tract column, `tract_len / period` times more often than
        // it can change, each time paying `period` emission lookups and a runtime `%` per unit
        // base. Here it costs `period²` lookups per row (≤ 36) once. The `k`-ascending
        // accumulation order is load-bearing: `f64` addition is not associative, so summing the
        // unit's bases in any other order would move the score's low bits.
        let mut unit_emit_by_phase = [UNREACHABLE; MAX_MOTIF_LEN];

        let UnitRobustScratch {
            rows,
            backpointers,
            column_plan,
        } = scratch;

        // Row 0 — no read base consumed. Every deletion here is a *leading* one (the read starts
        // later in the frame), so it takes the terminal cost regardless of the guard.
        {
            let row0 = &mut rows[0];
            row0[0] = [0.0, UNREACHABLE, UNREACHABLE, UNREACHABLE, UNREACHABLE];
            backpointers[0] = Backpointers::new([State::Match; STATES]);
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
                let (sd, sd_pred) = if column_plan[column].unit_deletable {
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
                backpointers[column] =
                    Backpointers::new([State::Match, State::Match, d_pred, State::Match, sd_pred]);
            }
        }

        for row_index in 1..=read_len {
            let read_base = bases[row_index - 1];
            let scores = self.emission.scores_for(read.quality_at(row_index - 1));
            let back_row = row_index * stride;

            let prev_slot = (row_index - 1) % ring_len;
            let slip_slot = row_index.checked_sub(period).map(|r| r % ring_len);
            let cur_slot = row_index % ring_len;

            // Resolve the whole-unit slip emission for this row, one entry per motif phase (see
            // `unit_emit_by_phase`). Only reachable once the row can be reached from `row - period`.
            if let Some(unit_start) = row_index.checked_sub(period) {
                for (phase, slot) in unit_emit_by_phase[..period].iter_mut().enumerate() {
                    let mut sum = 0.0;
                    // `motif_index` walks the motif from `phase`, wrapping once — both `phase` and
                    // `k` are below `period`, so their sum is below `2 * period` and one
                    // conditional subtraction is the whole modulo.
                    let mut motif_index = phase;
                    for k in 0..period {
                        let index = unit_start + k;
                        sum += self
                            .emission
                            .scores_for(read.quality_at(index))
                            .pick(bases[index], motif[motif_index]);
                        motif_index += 1;
                        if motif_index == period {
                            motif_index = 0;
                        }
                    }
                    *slot = sum;
                }
            }
            // A deletion on the last read row is a *trailing* one — the read has run out and the
            // rest of the frame must still be crossed. That is the `FromLeft` case, not the
            // boundary slide the guard is aimed at.
            let trailing = row_index == read_len;

            // Column 0 — a read base before any reference base: the leading-insertion route, how a
            // read that overhangs the frame's start is placed. Kept at the ordinary open cost for
            // the same reason the terminal deletions are.
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
            backpointers[back_row] = Backpointers::new([
                State::Match,
                ins0_pred,
                State::Match,
                State::Match,
                State::Match,
            ]);

            // The motif phase of the column about to be scored, `(column - left_flank_len) % period`
            // carried forward instead of divided out: it advances by one per column and wraps once
            // per unit, so the tract cells below index `unit_emit_by_phase` without a runtime `%`.
            // Seeded at `column == 1` (where `column - left_flank_len` may be negative, hence the
            // `+ period` before the reduction); only its values at tract columns are ever read.
            let mut phase = (1 + period - left_flank_len % period) % period;

            for column in 1..=reference_len {
                let emit = scores.pick(read_base, reference[column - 1]);
                let plan = column_plan[column];

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

                // Out-of-frame single-base insertion (a read base with no reference base). At the
                // frame's last column this is a *trailing* overhang — the read is longer than the
                // frame — so it takes the terminal cost, as the leading route at column 0 does.
                let ins_open = if column == reference_len {
                    plan.gap_open_terminal
                } else {
                    plan.gap_open
                };
                let (ins, ins_pred) = best_of(&[
                    (
                        ins_open + rows[prev_slot][column][State::Match.index()],
                        State::Match,
                    ),
                    (
                        self.costs.ln_gap_extend()
                            + rows[prev_slot][column][State::Insertion.index()],
                        State::Insertion,
                    ),
                ]);

                // Out-of-frame single-base deletion (a reference base with no read base). On the
                // last read row the terminal route applies, so a guarded open reopens.
                let del_open = if trailing {
                    plan.gap_open_terminal
                } else {
                    plan.gap_open
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
                // `i − period`, same column. Never guarded: real length variation is exactly what
                // this route is for.
                let (sins, sins_pred) = if let Some(slip_slot) = slip_slot {
                    if plan.in_tract {
                        let unit_emit = unit_emit_by_phase[phase];
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
                let (sdel, sdel_pred) = if plan.unit_deletable {
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
                    Backpointers::new([m_pred, ins_pred, del_pred, sins_pred, sdel_pred]);

                phase += 1;
                if phase == period {
                    phase = 0;
                }
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
            anchor: self.config.anchor,
        };
        Some(trace_back(
            final_state,
            backpointers,
            geometry,
            Aligned {
                read: bases,
                reference,
            },
        ))
    }
}

/// The matrix shape the traceback needs — bundled so the walk takes one `Copy` argument rather
/// than six positional `usize`s, several of them interchangeable (a transposition of
/// `read_len`/`reference_len` or of the two flank lengths would be a silent wrong measurement).
#[derive(Debug, Clone, Copy)]
struct MatrixGeometry {
    stride: usize,
    read_len: usize,
    reference_len: usize,
    left_flank_len: usize,
    right_flank_len: usize,
    period: usize,
    /// The evidence test the walk applies at the end.
    anchor: AnchorRule,
}

/// The two sequences the traceback compares, bundled so the anchor test can ask whether a `Match`
/// step really matched. Algorithm 4's walk never looks at the bases — its anchor test is purely
/// positional — so this is the one input this algorithm adds.
#[derive(Debug, Clone, Copy)]
struct Aligned<'a> {
    read: &'a [u8],
    reference: &'a [u8],
}

impl Aligned<'_> {
    /// Whether read base `read_index` and reference base `reference_index` are the *same* base —
    /// the evidence unit the anchor test counts.
    ///
    /// Case-insensitive, because neither a read nor a reference is guaranteed upper-case and a
    /// lower-case flank must not silently read as unanchored. An `N` on either side is never
    /// evidence: it is the absence of a base call, and counting it would let a run of `N`s anchor
    /// a flank the read never read.
    #[inline]
    fn bases_agree(&self, read_index: usize, reference_index: usize) -> bool {
        let read = self.read[read_index].to_ascii_uppercase();
        let reference = self.reference[reference_index].to_ascii_uppercase();
        read == reference && read != b'N'
    }
}

/// What one flank's traceback showed: how many read bases landed on it at all, and how many of
/// those, inside the junction window, actually agreed.
#[derive(Debug, Clone, Copy, Default)]
struct FlankEvidence {
    /// Read bases aligned against any column of this flank — agreeing or not.
    support: usize,
    /// Read bases inside the junction window whose base agreed with the reference.
    matches: usize,
}

impl FlankEvidence {
    /// Whether this flank is firm enough to anchor, under `rule`, at a flank of `flank_len` bases.
    /// Both thresholds are capped by the flank that exists, so a repeat at a contig edge is judged
    /// against the flank it has rather than made permanently unmeasurable.
    #[inline]
    fn is_firm(self, rule: AnchorRule, flank_len: usize) -> bool {
        self.matches >= rule.min_matches.min(flank_len)
            && self.support >= rule.min_support.min(flank_len)
    }
}

/// Turn algorithm 4's positional anchoring plus the two flanks' evidence into the final anchoring,
/// **demoting only ever downward and never past a lower bound**.
///
/// Three cases, and the third is the one that separates this from algorithm 5:
///
/// - Neither side crossed: nothing to decide, the read lies inside the repeat.
/// - Exactly one side crossed: the read is already a lower bound. Demoting the surviving side
///   would make it [`RepeatSpan::Unanchored`] — the read discarded — which is not a conservative
///   correction but a lost read. So a one-sided span is passed through untouched.
/// - Both sides crossed (algorithm 4 would say `Between`): each side must be firm. If neither is,
///   the span still keeps the better-evidenced side rather than collapsing to nothing.
#[inline]
fn resolve_anchors(
    crossed: (bool, bool),
    evidence: (FlankEvidence, FlankEvidence),
    rule: AnchorRule,
    flank_lens: (usize, usize),
) -> (bool, bool) {
    let (crossed_left, crossed_right) = crossed;
    if !(crossed_left && crossed_right) {
        return crossed;
    }
    let (left, right) = evidence;
    let firm_left = left.is_firm(rule, flank_lens.0);
    let firm_right = right.is_firm(rule, flank_lens.1);
    if firm_left || firm_right {
        return (firm_left, firm_right);
    }
    // Neither flank is firm, yet both were crossed. Keep the better-evidenced one so the read
    // survives as a lower bound; ties go left, matching the spec's determinism preference.
    if right.matches > left.matches
        || (right.matches == left.matches && right.support > left.support)
    {
        (false, true)
    } else {
        (true, false)
    }
}

/// Walk the traceback, read the tract off the two flank junctions, and **check that each flank
/// actually held it** — algorithm 4's walk plus the evidence counting the anchor test needs.
fn trace_back(
    final_state: State,
    backpointers: &[Backpointers],
    geometry: MatrixGeometry,
    aligned: Aligned<'_>,
) -> TractReadout {
    let MatrixGeometry {
        stride,
        read_len,
        reference_len,
        left_flank_len,
        right_flank_len,
        period,
        anchor,
    } = geometry;
    let left_junction = left_flank_len; // first tract reference base
    let right_junction = reference_len - right_flank_len; // first right-flank base
    let mut tract_start = 0usize;
    let mut tract_end = read_len;

    // The two junction-abutting windows, in reference columns: the last `window` bases of the left
    // flank, and the first `window` of the right.
    let left_window = left_junction.saturating_sub(anchor.window)..left_junction;
    let right_window = right_junction..(right_junction + anchor.window).min(reference_len);
    let mut left = FlankEvidence::default();
    let mut right = FlankEvidence::default();

    let mut i = read_len;
    let mut j = reference_len;
    let mut state = final_state;
    while i != 0 || j != 0 {
        let pred = backpointers[i * stride + j].get(state);
        match state {
            State::Match => {
                let consumed = j - 1;
                if consumed == left_junction {
                    tract_start = i - 1;
                }
                if right_flank_len > 0 && consumed == right_junction {
                    tract_end = i - 1;
                }
                // Support: a read base aligned against a flank column, agreeing or not.
                // Matches: the same, inside the junction window, and agreeing.
                if consumed < left_junction {
                    left.support += 1;
                    if left_window.contains(&consumed) && aligned.bases_agree(i - 1, consumed) {
                        left.matches += 1;
                    }
                } else if consumed >= right_junction {
                    right.support += 1;
                    if right_window.contains(&consumed) && aligned.bases_agree(i - 1, consumed) {
                        right.matches += 1;
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
                // A whole unit of read bases inserted — no reference consumed.
                i -= period;
            }
            State::SlipDeletion => {
                // A whole unit of reference bases deleted — no read consumed. The deleted 0-based
                // reference indices are `j - period ..= j - 1`. A contraction at the very start of
                // the tract deletes the first tract base, which **is** the left junction, so the
                // crossing must be recorded exactly as an ordinary deletion records it. The right
                // junction is the first right-flank base and lies beyond any tract-interior
                // deletion, so it cannot be consumed here; the check is kept defensively.
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

    // Algorithm 4's positional test — did the path cross the junction at all — then the evidence
    // test on top of it, which can only demote and never removes the last anchor.
    let crossed = (
        left_flank_len > 0 && tract_start != 0,
        right_flank_len > 0 && tract_end != read_len,
    );
    let (left_anchored, right_anchored) = resolve_anchors(
        crossed,
        (left, right),
        anchor,
        (left_flank_len, right_flank_len),
    );

    // A demoted junction leaves a *stale offset* behind, and the offset is load-bearing: a
    // `FromLeft` span must run to the end of the read, because that is what makes it a lower bound
    // on the allele. Keeping the crossing the anchor test just rejected would report a bound
    // shorter than the repeat the read actually showed — the same systematically short allele the
    // span type exists to prevent. So a rejected side falls back to the read's own edge, which is
    // what algorithm 4 would have left there had the path not crossed at all.
    if !left_anchored {
        tract_start = 0;
    }
    if !right_anchored {
        tract_end = read_len;
    }

    TractReadout {
        tract_start: tract_start as u64,
        tract_end: tract_end as u64,
        left_anchored,
        right_anchored,
    }
}

impl<E: Emission> BestPathAligner for SsrUnitRobustAligner<E> {
    type Scratch = UnitRobustScratch;
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

    /// Contraction-biased parameters — HipSTR's fitted values are, and the direction asymmetry is
    /// the whole point of the algorithm-4 family.
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
        let aligner = SsrUnitRobustAligner::new(PerQualityEmission::new());
        let context = RepeatContext {
            geometry,
            stutter: model,
        };
        let quality = vec![35u8; read.len()];
        let bases = ReadBases::try_new(read, &quality).expect("matched lengths");
        let mut scratch = UnitRobustScratch::new();
        aligner.align(bases, reference, context, &mut scratch)
    }

    /// **Mandatory property 1: a clean read of the reference measures the reference length.**
    #[test]
    fn a_clean_read_measures_the_reference_tract() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let span = measure(&reference, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
    }

    /// **The `chr3:33,877,690` shape — a substitution at a long homopolymer's first base,
    /// under the model the run actually uses.** The true allele spells the 11-base poly-A
    /// tract as `C` + 10 `A`s: same length, substituted first base. At 30× on HG002, 10 of
    /// 23 reads carry it and ng's candidate table held a bare 10-`A` run — the `C` pushed
    /// out of the tract and the length read one unit short (tract-accuracy program, L4).
    /// The flanks and tract are the reference's own bases at that locus; the model is
    /// `hipstr_shipped`, which is what `--defaults` hands this aligner.
    #[test]
    fn a_substitution_at_a_homopolymer_edge_stays_in_the_tract() {
        // The exact window the run hands the aligner: 15-base flanks (the default
        // bundle-threshold margin), the reference's own bases.
        let (reference, geometry) = frame(
            b"TTCTAGTTTTAGTTT",
            b"AAAAAAAAAAA",
            b"CTAAAAACCATTTTT",
            b"A",
        );
        let read = b"TTCTAGTTTTAGTTTCAAAAAAAAAACTAAAAACCATTTTT";
        let span = measure(read, &reference, &geometry, &StutterModel::hipstr_shipped());
        assert_eq!(span.measured_length(), Some(11), "span: {span:?}");
        let observed = span.observed_span().expect("a measured span");
        assert_eq!(
            &read[observed.start as usize..observed.end as usize],
            b"CAAAAAAAAAA"
        );
    }

    /// **The same locus, as the 30× reads actually spell it.** Every variant read at
    /// `chr3:33,877,690` carries the allele as a one-base deletion in the left flank's
    /// `TTT` run plus the `C`: `...TTCTAGTTTTAGTT | CAAAAAAAAAA | CTAAA...`. The honest
    /// account is a flank deletion, a mismatched `C` at the tract's first column, and ten
    /// matched `A`s — tract `CAAAAAAAAAA`, same length as the reference's. The rival path
    /// absorbs the `C` into the flank as a mismatch and prices the missing base as a
    /// whole-unit tract contraction, reporting a bare 10-`A` run — which is exactly what
    /// ng's candidate table held, and it costs the truth its candidacy at this locus.
    ///
    /// Ignored, deliberately red: the flank-side [`JunctionGuard`] makes the honest path's
    /// gap-open UNREACHABLE (the deletion sits 1–3 columns from the junction, inside the
    /// 7-column guard this flank gets), and even unguarded, a whole-unit slip (≈ −2.9 nats)
    /// under-prices the flank gap-open (≈ −10.4). The fix is the tract-accuracy program's
    /// L4 decision (`doc/devel/ng/research/tract_accuracy_program_report.md`); un-ignore
    /// with it.
    #[test]
    #[ignore = "L4: junction-adjacent real variation is inexpressible under the guard — fix pending the program's L4 decision"]
    fn a_flank_indel_beside_the_junction_does_not_eat_the_tract_edge() {
        let (reference, geometry) = frame(
            b"TTCTAGTTTTAGTTT",
            b"AAAAAAAAAAA",
            b"CTAAAAACCATTTTT",
            b"A",
        );
        let read = b"TTCTAGTTTTAGTTCAAAAAAAAAACTAAAAACCATTTTT";
        let span = measure(read, &reference, &geometry, &StutterModel::hipstr_shipped());
        assert_eq!(span.measured_length(), Some(11), "span: {span:?}");
        let observed = span.observed_span().expect("a measured span");
        assert_eq!(
            &read[observed.start as usize..observed.end as usize],
            b"CAAAAAAAAAA"
        );
    }

    /// A genuine whole-unit expansion is measured at its own length — the guard must not touch the
    /// slip route.
    #[test]
    fn a_whole_unit_expansion_is_measured() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let read = b"ACGTACGTCAGCAGCAGCAGCAGCAGTTGGTTGGAT"; // +2 units
        let span = measure(read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(18));
    }

    /// **Mandatory property 2: out-of-frame changes keep a route.** An interrupted repeat must be
    /// measured verbatim, not forced onto a tidy unit count — the guard is clamped precisely so
    /// this stays true.
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

    /// **Mandatory property 3 — the shared model, and direction asymmetry.**
    #[test]
    fn a_contraction_is_cheaper_to_open_than_an_expansion() {
        let slip = SlipCosts::from_model(&contraction_biased());
        assert!(
            slip.open_contraction > slip.open_expansion,
            "a contraction-biased model must make contraction the cheaper slip to open"
        );
    }

    /// The slip costs reconstruct the stutter model's own probabilities exactly.
    #[test]
    fn the_slip_costs_reconstruct_the_stutter_probability() {
        let model = contraction_biased();
        let slip = SlipCosts::from_model(&model);
        let period = std::num::NonZeroU8::new(3).unwrap();
        let ln_same_length_share = model.same_length_share().ln();

        for n in 1..=5i64 {
            let reconstructed = slip.open_expansion + (n - 1) as f64 * slip.extend;
            let expected = model.probability(n * 3, period).ln() - ln_same_length_share;
            assert!((reconstructed - expected).abs() < 1e-12);
            let reconstructed = slip.open_contraction + (n - 1) as f64 * slip.extend;
            let expected = model.probability(-n * 3, period).ln() - ln_same_length_share;
            assert!((reconstructed - expected).abs() < 1e-12);
        }
    }

    /// **The reason the anchor test exists.** A long allele that outruns the read leaves no right
    /// flank in the read, yet the frame still has one to account for, and spending the read's last
    /// repeat bases as mismatched "matches" against it is cheaper than deleting it. Algorithm 4
    /// reads that crossing as an anchor and reports a fabricated complete length.
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
        let observed = span
            .observed_span()
            .expect("a lower bound still has a span");
        assert_eq!(observed.end, read.len() as u64);
    }

    /// The mirror case: a read that begins inside a long tract anchors only on the right.
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

    /// **A real 1 bp indel in the body of a flank must not derail the measurement.** This is what
    /// a flat flank ban gives up, and the whole reason the guard is a window: the tract length is
    /// unchanged, so the correct answer is the reference tract length.
    #[test]
    fn an_indel_in_the_body_of_a_flank_leaves_the_tract_measured() {
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
                    Some(ref_tract.len() as u64),
                    "flank indel (right_side={right_side}, insert={insert}) moved the tract"
                );
            }
        }
    }

    /// **The demotion must not fire on real reads.** A genuinely spanning read with a miscall in
    /// the junction window itself still measures its allele — the test is a quorum, not an exact
    /// run.
    #[test]
    fn a_miscall_beside_each_junction_does_not_demote_a_spanning_read() {
        let (reference, geometry) = frame(b"ACGTACGT", b"CAGCAGCAGCAG", b"TTGGTTGGAT", b"CAG");
        let mut read = reference.clone();
        read[7] = b'A'; // last base of the left flank, miscalled
        read[20] = b'A'; // first base of the right flank, miscalled
        let span = measure(&read, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
    }

    /// **The rule that separates this from algorithm 5: a partial is never demoted to nothing.**
    /// A one-sided crossing is passed through untouched however thin its evidence, and a
    /// two-sided one always keeps at least one side.
    #[test]
    fn a_partial_is_never_demoted_to_unanchored() {
        let rule = AnchorRule::default();
        let thin = FlankEvidence {
            support: 0,
            matches: 0,
        };
        // One-sided crossings pass through, evidence or not.
        assert_eq!(
            resolve_anchors((true, false), (thin, thin), rule, (25, 25)),
            (true, false)
        );
        assert_eq!(
            resolve_anchors((false, true), (thin, thin), rule, (25, 25)),
            (false, true)
        );
        // Neither crossed: nothing to decide.
        assert_eq!(
            resolve_anchors((false, false), (thin, thin), rule, (25, 25)),
            (false, false)
        );
        // Both crossed but neither firm: the better-evidenced side survives, never nothing.
        let better = FlankEvidence {
            support: 4,
            matches: 1,
        };
        assert_eq!(
            resolve_anchors((true, true), (thin, better), rule, (25, 25)),
            (false, true)
        );
        assert_eq!(
            resolve_anchors((true, true), (better, thin), rule, (25, 25)),
            (true, false)
        );
        // ...and a tie goes left, deterministically.
        assert_eq!(
            resolve_anchors((true, true), (thin, thin), rule, (25, 25)),
            (true, false)
        );
    }

    /// A two-sided crossing with firm evidence on both sides stays a measurement.
    #[test]
    fn a_firm_two_sided_crossing_stays_a_measurement() {
        let rule = AnchorRule::default();
        let firm = FlankEvidence {
            support: 20,
            matches: 5,
        };
        assert_eq!(
            resolve_anchors((true, true), (firm, firm), rule, (25, 25)),
            (true, true)
        );
    }

    /// A locus whose flanks are shorter than the quorum is judged against the flank it has —
    /// otherwise every contig-edge locus would become permanently unmeasurable.
    #[test]
    fn a_flank_shorter_than_the_quorum_can_still_anchor() {
        let (reference, geometry) = frame(b"T", b"CAGCAGCAGCAG", b"T", b"CAG");
        let span = measure(&reference, &reference, &geometry, &contraction_biased());
        assert_eq!(span.measured_length(), Some(12));
    }

    /// **The guard is clamped so the tract keeps an interior**, and so a flank keeps a body. A
    /// short tract must not have its whole out-of-frame route closed.
    #[test]
    fn the_guard_is_clamped_to_leave_an_interior_and_a_flank_body() {
        let guard = JunctionGuard::default();
        // A twelve-base, four-unit period-3 tract: (12 − 3)/2 = 4 columns a side, not six.
        assert_eq!(guard.tract_width(12, 3), 4);
        // A tract no longer than one unit has no interior to guard at all.
        assert_eq!(guard.tract_width(3, 3), 0);
        assert_eq!(guard.tract_width(0, 3), 0);
        // A long tract takes the configured width.
        assert_eq!(guard.tract_width(60, 3), guard.tract_columns);
        // Flanks: never more than half, so a real indel in the body always keeps a route.
        assert_eq!(guard.flank_width(25), guard.flank_columns);
        assert_eq!(guard.flank_width(8), 4);
        assert_eq!(guard.flank_width(1), 0);
        assert_eq!(guard.flank_width(0), 0);
    }

    /// The quorum is counted inside the window, so it can only be met if it fits there — and a
    /// quorum of zero is algorithm 4's test with extra steps.
    #[test]
    fn the_anchor_quorum_fits_inside_the_anchor_window() {
        let rule = AnchorRule::default();
        assert!(rule.min_matches <= rule.window);
        assert!(rule.min_matches > 0);
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

    /// Scratch reuse across a size drop must not leak — the ring is grown, never cleared.
    #[test]
    fn scratch_reuse_does_not_leak_across_periods() {
        let aligner = SsrUnitRobustAligner::new(PerQualityEmission::new());
        let model = contraction_biased();

        let loci: &[(&[u8], &[u8])] = &[
            (b"GTTGTG", b"GTTGTGGTTGTGGTTGTGGTTGTGGTTGTG"),
            (b"A", b"AAAA"),
            (b"CAG", b"CAGCAGCAGCAG"),
            (b"CA", b"CACACACACACA"),
        ];

        let mut reused_scratch = UnitRobustScratch::new();
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
            let fresh = aligner.align(bases, &reference, context, &mut UnitRobustScratch::new());
            assert_eq!(
                reused,
                fresh,
                "scratch reuse leaked across periods at motif {:?}",
                std::str::from_utf8(motif).unwrap()
            );
        }
    }

    /// Multiple periods work — a period-1 (homopolymer) and a period-2 locus, since the slip
    /// stride is the period and period 1 is where the indel deficit lives.
    #[test]
    fn homopolymer_and_dinucleotide_loci_are_measured() {
        let (reference, geometry) = frame(b"ACGTACGT", b"AAAAAA", b"TTGGTTGGAT", b"A");
        assert_eq!(
            measure(&reference, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(6)
        );
        let read = b"ACGTACGTAAAAAAATTGGTTGGAT";
        assert_eq!(
            measure(read, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(7)
        );

        let (reference, geometry) = frame(b"ACGTACGT", b"CACACACA", b"TTGGTTGGAT", b"CA");
        assert_eq!(
            measure(&reference, &reference, &geometry, &contraction_biased()).measured_length(),
            Some(8)
        );
    }
}
