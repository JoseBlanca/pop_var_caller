//! The fold's witnessed set — what one read saw of one *record*, in reference coordinates.
//!
//! A file of its own rather than a corner of `open_record.rs`, and for the reason
//! [`witness.rs`](crate::ng::locus_generation) exists: a private field is private **to a
//! module**, and `open_record.rs` is where `apply_events_into`, `fold_read_into_record`,
//! `finalise` and `witness_of` live — so a field private *there* is reachable from exactly
//! the code most tempted to bypass it. Here, the only way to build one of these is through
//! a constructor that canonicalises (Milestone B review, owner 2026-07-30).
//!
//! Visibility is unchanged by the move: `pub(super)`, so the type stays the walk's own and
//! nothing outside `pileup` can name it (arch §2).

use smallvec::SmallVec;

use crate::ng::locus_generation::witness::canonicalise_runs;

/// The reference positions one read witnessed inside a record — the fold-time counterpart
/// of [`WitnessedLocusPositions`](crate::ng::locus_generation::WitnessedLocusPositions), and
/// what C1 put in place of `RefSpan` — the one-run, inclusive-ended pair the fold used to
/// return, deleted in the same commit.
///
/// Canonical half-open runs in **reference** coordinates, sorted, non-empty, and separated
/// by at least one position the read did not witness. It carries `RefSpan`'s invariant
/// verbatim: **there is no empty set** — a read that witnessed nothing inside a record does
/// not fold into it at all, so "no positions" is the absence of a `WitnessedRefPositions`,
/// never a degenerate one.
///
/// # Two types, one per axis — and why this one is not in `witness.rs`
///
/// The fold works in reference coordinates and `finalise` resolves into locus coordinates
/// against the record's **final** footprint. Those are two axes, and they get two types for
/// the reason that minted `LocusLen`: same-shaped quantities that can be transposed get two
/// types so the compiler refuses the mix (arch §3). This one is the *fold's*, `pub(super)`
/// and private to the walk, so it lives with the walk; only the locus-axis vocabulary is
/// public (arch *Module home*). The canonical form itself is shared — one
/// [`canonicalise_runs`](crate::ng::locus_generation::witness::canonicalise_runs), so the
/// two types cannot drift on what "canonical" means.
///
/// # Half-open, where `RefSpan` was inclusive
///
/// `RefSpan` is inclusive of both ends, matching `alignment_start`/`alignment_end`. A *set*
/// is held half-open instead, because merging is where an off-by-one silently changes the
/// answer: with half-open runs, "touching" is `start == end` and needs no `+ 1` anywhere,
/// and the fold already accumulates half-open before converting on the way out
/// (`apply_events_into`'s `witnessed` is `(start, end)` exclusive today). `witness_of` is
/// the only conversion off this type, so the convention stays inside the walk.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct WitnessedRefPositions(SmallVec<[(u32, u32); 2]>);

/// The scratch a fold accumulates one read's runs into, so the buffer is the caller's and
/// the per-read cost is a swap rather than an allocation (arch §2).
pub(super) type WitnessedRefRuns = SmallVec<[(u32, u32); 2]>;

impl WitnessedRefPositions {
    /// From runs in **any** order, normalised: sorted by start, then touching and
    /// overlapping runs merged. Half-open `[start, end)`. `None` for an empty iterator or
    /// any run with `start >= end`.
    ///
    /// The fold does not use this — it owns a buffer and goes through
    /// [`take_from`](Self::take_from) / [`refill_from`](Self::refill_from) so a re-fold
    /// costs a swap rather than an allocation. This is the shape a *test* wants, where the
    /// runs are written out literally, and the `expect` is what makes a caller who reaches
    /// for it on the hot path notice.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the fold owns a buffer and goes through take_from/refill_from"
        )
    )]
    pub(super) fn from_half_open_runs(runs: impl IntoIterator<Item = (u32, u32)>) -> Option<Self> {
        canonicalise_runs(runs.into_iter().collect()).map(Self)
    }

    /// Take a caller-owned scratch buffer's runs, canonicalise them, and leave the buffer
    /// empty for the next read. `None` — and the buffer untouched — when the read witnessed
    /// nothing, which is the case the fold answers with no observation at all.
    ///
    /// This is the shape arch §2 asks for: "the runs are accumulated into a buffer the
    /// caller owns and the callee clears, like `allele_seq`, so a fold allocates nothing
    /// per read". Its sibling [`refill_from`](Self::refill_from) is what keeps that true
    /// across *re*-folds.
    pub(super) fn take_from(buf: &mut WitnessedRefRuns) -> Option<Self> {
        canonicalise_runs(std::mem::take(buf)).map(Self)
    }

    /// Replace this set from a caller-owned scratch buffer, **swapping** rather than
    /// moving, so the buffer inherits this set's storage and a fold that has spilled once
    /// allocates nothing thereafter. `false` — with both this set and the buffer unchanged
    /// — when the runs describe nothing.
    ///
    /// The fold needs this because a read's witness is rebuilt on **every widen**, not once
    /// per read: `refold_live_reads` re-places every live folded read against the wider
    /// window. A witness of three or more runs owns a heap buffer, so a move-based refill
    /// would allocate per (read × widen) on exactly the multi-junction RNA-seq case this
    /// milestone exists for (Milestone B review).
    pub(super) fn refill_from(&mut self, buf: &mut WitnessedRefRuns) -> bool {
        let Some(canonical) = canonicalise_runs(std::mem::take(buf)) else {
            // `take` emptied the buffer; give it back, so "nothing witnessed" leaves both
            // sides exactly as they were.
            return false;
        };
        *buf = std::mem::replace(&mut self.0, canonical);
        buf.clear();
        true
    }

    /// The runs, in canonical order, half-open `[start, end)` in reference coordinates.
    ///
    /// **The only way to read the set**, deliberately. It used to offer `start()` and
    /// `end_exclusive()` as well — the old `RefSpan`'s two fields — and the Milestone B
    /// review showed what that invites: `witness_of` computes `past_last` as
    /// `witnessed.end.saturating_add(1)` against the *inclusive* `RefSpan`, and the cheapest
    /// C1 port substitutes `end_exclusive()` there. It type-checks, adds a position to every
    /// witness, and reinstates the span-swallows-the-hole behaviour this milestone exists to
    /// remove — silently, because `positions_covered` is a third of an observation's
    /// identity. Both accessors were dead code, so they are gone; C1 walks the runs.
    pub(super) fn runs(&self) -> impl ExactSizeIterator<Item = (u32, u32)> + '_ {
        self.0.iter().copied()
    }

    /// How many reference positions the read witnessed — **the positions, not the span**.
    /// The number under the obvious name, at the point where the wrong one (last end minus
    /// first start) is tempting.
    ///
    /// The fold has no use for it: `witness_of` resolves the runs against the *footprint*,
    /// which is a different number from the runs' own total, and nothing else on this axis
    /// wants a count. It stays because the tests below are what pin the merge against the
    /// enclosing span, and it is the reference-axis half of the number the locus axis emits.
    #[cfg_attr(
        not(test),
        expect(
            dead_code,
            reason = "the fold measures against the footprint, not the runs' own total"
        )
    )]
    pub(super) fn positions_covered(&self) -> u64 {
        self.0
            .iter()
            .map(|(start, end)| u64::from(*end) - u64::from(*start))
            .sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // -----------------------------------------------------------------
    // B3 — the fold's witnessed set, in reference coordinates.
    // -----------------------------------------------------------------

    /// **`RefSpan`'s invariant, carried over verbatim: there is no empty set.** A read that
    /// witnessed nothing inside a record does not fold into it, so "no positions" is the
    /// absence of a `WitnessedRefPositions` and never a degenerate one — the same statement
    /// `RefSpan`'s doc makes, on a type that now holds several runs.
    #[test]
    fn a_witnessed_ref_set_is_never_empty_and_never_holds_an_empty_run() {
        assert!(WitnessedRefPositions::from_half_open_runs([]).is_none());
        assert!(
            WitnessedRefPositions::from_half_open_runs([(20, 20)]).is_none(),
            "half-open, so start == end covers nothing",
        );
        assert!(
            WitnessedRefPositions::from_half_open_runs([(5, 9), (20, 20)]).is_none(),
            "an empty run among good ones is rejected rather than dropped",
        );
    }

    /// The reference axis normalises exactly as the locus axis does — same rule, one
    /// implementation, asserted here so the *reference* type is not taken on trust.
    #[test]
    fn witnessed_ref_runs_sort_and_merge_and_keep_their_holes() {
        let merged = WitnessedRefPositions::from_half_open_runs([(20, 24), (5, 9), (9, 12)])
            .expect("non-empty");
        assert_eq!(
            merged.runs().collect::<Vec<_>>(),
            vec![(5, 12), (20, 24)],
            "out of order, and the two that touch at 9 are one run",
        );
        assert_eq!(
            merged.positions_covered(),
            11,
            "7 positions in the merged run and 4 in the second — not the 19 the enclosing \
             span would give",
        );
    }

    /// **Containment**, the case that separates "merge" from "extend": a run wholly inside
    /// another adds no positions and must not survive as a second run, and the merged end
    /// must stay the *outer* one. The locus axis has a fixture for it; the reference axis
    /// needs its own, because a normaliser that took the inner end instead would lose
    /// witnessed positions while leaving every shape property — sorted, disjoint, non-empty
    /// — intact (Milestone B review).
    #[test]
    fn a_run_contained_in_another_neither_survives_nor_shortens_it() {
        let contained =
            WitnessedRefPositions::from_half_open_runs([(5, 24), (9, 12)]).expect("non-empty");
        assert_eq!(contained.runs().collect::<Vec<_>>(), vec![(5, 24)]);
        assert_eq!(contained.positions_covered(), 19);
    }

    /// The buffer-reusing constructors arch §2 asks for: the fold accumulates into a buffer
    /// it owns, and the callee leaves it empty. `refill_from` **swaps**, so the buffer
    /// inherits the old set's storage rather than dropping it — which is what makes a
    /// re-fold free after the first spill.
    #[test]
    fn taking_and_refilling_leave_the_callers_buffer_empty_and_reusable() {
        let mut buf: WitnessedRefRuns = SmallVec::from_slice(&[(9, 12), (5, 9)]);
        let mut set = WitnessedRefPositions::take_from(&mut buf).expect("two runs, one merged");
        assert_eq!(set.runs().collect::<Vec<_>>(), vec![(5, 12)]);
        assert!(buf.is_empty(), "the callee clears the caller's buffer");

        buf.extend_from_slice(&[(40, 44), (20, 24)]);
        assert!(set.refill_from(&mut buf));
        assert_eq!(set.runs().collect::<Vec<_>>(), vec![(20, 24), (40, 44)]);
        assert!(buf.is_empty(), "and again on a refill");

        // Nothing witnessed: both sides unchanged, and the caller's buffer still usable.
        buf.push((7, 7));
        assert!(!set.refill_from(&mut buf), "an empty run is not a set");
        assert_eq!(
            set.runs().collect::<Vec<_>>(),
            vec![(20, 24), (40, 44)],
            "a failed refill must not leave a half-written set behind",
        );
        assert!(WitnessedRefPositions::take_from(&mut SmallVec::new()).is_none());
    }

    /// **A hole is what the type exists to hold**, and the reference axis is where the fold
    /// first sees one: a spliced read's two exons, or an interior `N`. One run per exon,
    /// not one span swallowing the intron.
    #[test]
    fn two_runs_with_a_hole_between_them_stay_two_runs() {
        let spliced = WitnessedRefPositions::from_half_open_runs([(100, 103), (118, 121)])
            .expect("non-empty");
        assert_eq!(spliced.runs().len(), 2);
        assert_eq!(
            spliced.runs().collect::<Vec<_>>(),
            vec![(100, 103), (118, 121)]
        );
        assert_eq!(
            spliced.positions_covered(),
            6,
            "6 positions witnessed across a 21-position span — which is why the fold \
             cannot keep answering with a span",
        );
    }
}
