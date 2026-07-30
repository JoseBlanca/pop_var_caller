//! What one read saw of one locus — the witness vocabulary.
//!
//! [`ReadWitness`] is the answer to "how much of this locus did this read witness",
//! [`WitnessedLocusPositions`] is the set of positions that will replace its single run at
//! Milestone C, and [`LocusLen`] is the axis both are measured on. [`canonicalise_runs`] is
//! the one normaliser every witnessed-set type shares — including the fold's reference-axis
//! `WitnessedRefPositions`, which lives with the walk.
//!
//! They live in their own file rather than in `mod.rs` because
//! `WitnessedLocusPositions`'s invariant — one set of witnessed positions has exactly one
//! representation, so two reads that witnessed the same thing share one observation —
//! fails *silently* when it fails, and a private field needs one module boundary to be
//! private to (arch *Module home*).
//!
//! The three types are re-exported from [`locus_generation`](super), so no *external*
//! consumer's import path names this module; `canonicalise_runs` is crate-internal and is
//! imported by path, by `pileup::open_record`.

use smallvec::{Array, SmallVec};

use crate::ng::types::GenomeRegion;

/// Sort and merge `runs` into the canonical form every witnessed-set type holds: sorted by
/// start, non-empty, and separated by at least one position. `None` for an empty input or
/// any empty run, which is the "no empty span" invariant every one of them carries.
///
/// **The buffer flows through the return value** so that ignoring the answer cannot
/// compile. It reported success as a `bool` and canonicalised in place until the Milestone
/// B review pointed out what that admits: `canonicalise_runs(&mut runs); Some(Self(runs))`
/// builds a degenerate set — empty, or holding empty runs — and no lint says a word.
/// Taking the buffer by value and handing it back closes that off, and the callers that
/// reuse a buffer across reads swap it back (see `WitnessedRefPositions::refill_from`).
///
/// **One normaliser, several types.** The locus axis and the reference axis are
/// deliberately separate types so the compiler refuses to mix them (arch §3), but the
/// *canonical form* is one rule, and two copies of it could drift while both still
/// compiled — which is exactly what happened to `witness_order` and cost Milestone A a
/// review finding. The types stay distinct; the rule has one home.
pub(in crate::ng::locus_generation) fn canonicalise_runs<T, A>(
    mut runs: SmallVec<A>,
) -> Option<SmallVec<A>>
where
    T: Copy + Ord,
    A: Array<Item = (T, T)>,
{
    if runs.is_empty() || runs.iter().any(|(start, end)| start >= end) {
        return None;
    }
    runs.sort_unstable();
    // Merge left to right, writing back into the same buffer. `start <= open_end` covers
    // both cases the invariant forbids: overlapping (`start < open_end`) and merely
    // touching (`start == open_end`), which are one run of witnessed positions and must
    // have one spelling.
    let mut kept = 0usize;
    for index in 1..runs.len() {
        let (start, end) = runs[index];
        let open_end = runs[kept].1;
        if start <= open_end {
            runs[kept].1 = open_end.max(end);
        } else {
            kept += 1;
            runs[kept] = (start, end);
        }
    }
    runs.truncate(kept + 1);
    Some(runs)
}

/// The locus positions one read witnessed — **a set, in locus coordinates**, and never
/// empty (a read that witnessed nothing inside a footprint does not fold into it).
///
/// Held as **canonical half-open runs**: sorted by start, each non-empty, and separated by
/// at least one position the read did not witness. Two runs that touch or overlap are one
/// run, always.
///
/// # The invariant is the whole point, and it fails silently
///
/// An observation's identity is `(bases, read_witness, read_group)`, and observations are
/// merged by comparing it. **If two equal sets could have two representations, two reads
/// with the same witness would stop sharing an observation** — not as an error, but as
/// extra observations, i.e. a wrong number with nothing to catch it (spec §3.3). So the
/// field is private, every constructor normalises, and `Eq`/`Hash` are derived over a
/// representation that admits exactly one spelling per set. `SmallVec` compares and hashes
/// as a slice, so a set that fits inline and the same set after a spill are equal.
///
/// # Two runs inline
///
/// Every witness in 225 million DNA-seq event-folds is one run, and the RNA-seq case this
/// change exists for is two (spec §8), so the common path allocates nothing. A bitmask was
/// rejected on the same evidence: it would pay 16 bytes to say "positions 0..40" for a case
/// that is one run in every observation measured (arch §4). The encoding is private behind
/// [`runs`](Self::runs) so it can still move.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct WitnessedLocusPositions(SmallVec<[(u16, u16); 2]>);

impl WitnessedLocusPositions {
    /// From runs in **any** order, normalised: sorted by start, then touching and
    /// overlapping runs merged.
    ///
    /// Runs are half-open `[start, end)` — which is what the name says, and it says it
    /// because the type's other constructor speaks the other convention. `None` when there
    /// is nothing to describe: an empty iterator, or any run with `start >= end`, since a
    /// run of no positions is not something a read can witness and silently dropping it
    /// would let a caller build a set it did not mean.
    pub fn from_half_open_runs(runs: impl IntoIterator<Item = (u16, u16)>) -> Option<Self> {
        canonicalise_runs(runs.into_iter().collect()).map(Self)
    }

    /// The one-run case — what the STR path mints and the common generic one — taking the
    /// run as **the offset it starts at and how many positions it covers**.
    ///
    /// # Why the name carries the convention
    ///
    /// Both spellings are a pair of `u16` locus positions, and the misread is not rejected:
    /// handing `from_half_open_runs` an offset/length pair produces a perfectly canonical
    /// set that is silently the *wrong* set — `(4, 5)` covers one position where the caller
    /// meant five. Only the reversed spelling errors. That is the hazard [`LocusLen`]'s own
    /// doc says it was minted to remove ("two `u16`s in a row with no way to tell them
    /// apart"), so the two constructors say which language they speak at every call site
    /// (Milestone B review). Newtypes for the two components would close it completely, and
    /// are the right move when `ReadWitness::Partial`'s fields are revisited — they are the
    /// same two `u16`s, and sealing that variant is deferred (spec §6).
    ///
    /// `None` if the run covers no positions, or if it would end past `u16::MAX` — the
    /// out-of-range case the spec asks to fail rather than clamp (§3.4), since a clamp here
    /// would silently shorten a witness.
    pub fn one_run_from_offset_and_length(
        offset_in_locus: u16,
        positions_covered: u16,
    ) -> Option<Self> {
        let end = offset_in_locus.checked_add(positions_covered)?;
        Self::from_half_open_runs([(offset_in_locus, end)])
    }

    /// The runs, in canonical order, half-open `[start, end)` in locus coordinates.
    pub fn runs(&self) -> impl ExactSizeIterator<Item = (u16, u16)> + '_ {
        self.0.iter().copied()
    }

    /// How many locus positions in total — what `num_obs_along_locus` sums over.
    ///
    /// **The positions, not the span.** Two runs with a gap between them cover fewer
    /// positions than the distance from the first start to the last end, and the difference
    /// is exactly what this type exists to record. `u32` because the runs can together
    /// exceed a `u16` only if they overlap, which they never do — but the sum is taken
    /// without that assumption.
    pub fn positions_covered(&self) -> u32 {
        self.0
            .iter()
            .map(|(start, end)| u32::from(*end) - u32::from(*start))
            .sum()
    }

    /// Whether the witness starts at the locus's left border — a **prefix** constraint on
    /// the allele, unchanged in meaning from the one-run form.
    pub fn is_flush_left(&self) -> bool {
        self.0.first().is_some_and(|(start, _)| *start == 0)
    }

    /// Whether the witness reaches the locus's right border — a **suffix** constraint.
    ///
    /// `>=` rather than `==` for the same reason the one-run predicate uses it: a producer's
    /// reach is measured in read bases, which diverge from locus positions under stutter, so
    /// a run may be handed a length the locus does not have.
    pub fn is_flush_right(&self, locus_len: LocusLen) -> bool {
        self.0
            .last()
            .is_some_and(|(_, end)| *end >= locus_len.get())
    }
}

/// How much of a locus a single read spanned — **one read's span, not depth**.
///
/// `Complete` means the read reached **both** borders of the locus; anything else is
/// the one **run** of locus positions the read actually witnessed, in **locus**
/// coordinates (spec §3). A partial run is a *censored* observation: the sequence is
/// at least this long, but not how long.
///
/// **One `Partial` run replaces the earlier `PartialLeft`/`PartialRight` pair
/// (owner, 2026-07-28).** Two side-tagged variants cannot describe what a read
/// witnesses once the *events*, not the alignment span, define it: a read can be blind
/// in the **middle** of a footprint (an interior `N`, a ref-skip) or blind at either
/// end, and a widened record can be wider than the read on both sides. Prefix-versus-
/// suffix survives as a derivation — [`is_flush_left`](Self::is_flush_left) /
/// [`is_flush_right`](Self::is_flush_right) — so the STR path's "a prefix and a suffix
/// are different constraints" is preserved, not lost. `Complete` is kept rather than
/// folded in: it is the overwhelmingly common case and it keeps
/// [`complete_observations`](super::SampleLocusObservations::complete_observations) a cheap
/// equality instead of arithmetic against the footprint.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ReadWitness {
    /// The read reached both borders of the locus.
    Complete,
    /// The stretch the read did witness, in **locus positions** — the axis `bases` is
    /// not on (that is allele content, in read coordinates).
    ///
    /// **Named `Partial`, not `Observed`.** Next to `Complete`, "observed" is not a
    /// contrast — a complete witness was observed too — and once the enum itself says
    /// *witness*, the word adds nothing. `Partial` says the one thing that separates
    /// the two (spec §3.1).
    ///
    /// **The fields are public, and the clamping in [`from_left`](Self::from_left) /
    /// [`from_right`](Self::from_right) is therefore a convention rather than a type
    /// invariant.** Left that way deliberately (2026-07-28): private fields would prove
    /// only that a run had been clamped against *some* [`LocusLen`], and nothing ties
    /// that length to the locus the run is finally attached to — `ReadWitness` cannot
    /// know its own locus. So the real check has to live where the region is in hand,
    /// which is `num_obs_along_locus`, and it does.
    ///
    /// Revisit when the **generic** path mints its first run: it needs runs flush with
    /// neither border (a read blind in the middle of a footprint), which neither
    /// constructor expresses, so the full constructor set — and with it the case for
    /// sealing the variant — is only knowable then. Building it now would be designing
    /// against one producer and guessing at the second.
    Partial {
        /// Locus positions between the locus's left border and the first one
        /// witnessed. `0` = flush with the left border, i.e. a prefix constraint.
        offset_in_locus: u16,
        /// How many locus positions were witnessed, running from `offset_in_locus`.
        positions_covered: u16,
    },
}

/// A locus's length in reference positions — the axis a [`ReadWitness`] run lives on.
///
/// **A newtype because the alternative is silently wrong.** Every constructor and
/// predicate on `ReadWitness` takes a covered extent *and* a locus length, both counts
/// of locus positions and both formerly `u16` — so `from_left(10, 4)` and
/// `from_left(4, 10)` each compiled, and the clamping the constructors do for their own
/// good would then hide the transposition rather than surface it. Two `u16`s in a row
/// with no way to tell them apart is exactly the shape that produces a wrong depth and
/// no panic.
///
/// It also gives the saturating cast **one** home. The count arrives as a `u64` (a
/// region length, a tract length) and has to be narrowed; doing that at each call site
/// spread the same `.min(u16::MAX as u64) as u16` across the mint and every dump tool.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocusLen(u16);

impl LocusLen {
    /// From a count of reference positions, saturating at `u16::MAX`.
    ///
    /// Saturation rather than an error: a locus longer than 65,535 positions cannot be
    /// described by a run either, and the clamp keeps the derivation total. A tract that
    /// long is a satellite, which the caller does not handle.
    pub fn from_positions(positions: u64) -> Self {
        Self(positions.min(u16::MAX as u64) as u16)
    }

    /// The length of `region` — the canonical source once a locus exists, and what
    /// [`SampleLocusObservations::locus_len`](super::SampleLocusObservations::locus_len) returns.
    pub fn of_region(region: GenomeRegion) -> Self {
        Self::from_positions(region.len())
    }

    pub fn get(self) -> u16 {
        self.0
    }
}

impl ReadWitness {
    /// A run flush with the locus's **left** border, `positions_covered` long — the
    /// old `PartialLeft(n)`.
    ///
    /// `locus_len` clamps the run into the locus. It is needed because a producer's
    /// reach can be measured in *read* bases, which diverge from locus positions under
    /// stutter; a run must never claim positions the locus does not have.
    pub fn from_left(positions_covered: u16, locus_len: LocusLen) -> Self {
        Self::Partial {
            offset_in_locus: 0,
            positions_covered: positions_covered.min(locus_len.get()),
        }
    }

    /// A run flush with the locus's **right** border, `positions_covered` long — the
    /// old `PartialRight(n)`.
    ///
    /// Clamped first, then the offset derived, so the subtraction cannot underflow:
    /// computing `locus_len - positions_covered` on an over-long reach would wrap to a
    /// huge offset, or — with a saturating subtraction — silently relabel a
    /// right-anchored read as a left-anchored one.
    ///
    /// **One consequence of the encoding, and it reaches further than labelling.** A run
    /// that covers the whole locus is flush with *both* borders, so once
    /// `positions_covered >= locus_len` this returns exactly what
    /// [`from_left`](Self::from_left) would. Three things follow, and the third is a
    /// behaviour change:
    ///
    /// 1. the run is reported as flush left by [`is_flush_left`](Self::is_flush_left);
    /// 2. every label derived from flushness calls it a *left* partial;
    /// 3. **it shares a bucket key with a left-flush run of the same bases, so the two
    ///    merge into one observation** — where `PartialLeft(n)` and `PartialRight(n)` kept them
    ///    apart.
    ///
    /// That is reachable on the STR path, where the reach is measured in *read* bases:
    /// an allele longer than the reference tract gives a reach past the locus length.
    /// It is arguably the right answer — identical constraints are one observation, and a read
    /// that witnessed every position is constrained from neither side — but it is not
    /// the pre-reshape answer, and the plan's stated equivalence
    /// `PartialRight(n) ⇔ Partial { len - n, n }` stops holding at `n = len`. Pinned by
    /// `ssr::tally::tests::an_expanded_allele_merges_the_two_sides_into_one_observation`.
    pub fn from_right(positions_covered: u16, locus_len: LocusLen) -> Self {
        let covered = positions_covered.min(locus_len.get());
        Self::Partial {
            offset_in_locus: locus_len.get() - covered,
            positions_covered: covered,
        }
    }

    /// Whether the run starts at the locus's left border — a **prefix** constraint on
    /// the allele. Always true of `Complete`.
    pub fn is_flush_left(&self) -> bool {
        match self {
            Self::Complete => true,
            Self::Partial {
                offset_in_locus, ..
            } => *offset_in_locus == 0,
        }
    }

    /// Whether the run ends at the locus's right border — a **suffix** constraint.
    /// Always true of `Complete`.
    pub fn is_flush_right(&self, locus_len: LocusLen) -> bool {
        match self {
            Self::Complete => true,
            Self::Partial {
                offset_in_locus,
                positions_covered,
            } => offset_in_locus.saturating_add(*positions_covered) >= locus_len.get(),
        }
    }

    /// A **total order over witnesses**: complete first, then partial runs by where they
    /// sit in the locus and then by how far they run. Both generators sort their
    /// observations by it, so two loci with the same evidence present it the same way.
    ///
    /// **A named key rather than an `Ord` impl**, so the convention is opt-in at the sort
    /// site and a witness is not silently comparable everywhere. That much is unchanged;
    /// what changed is *where it lives*. `pileup/open_record.rs` and `ssr.rs`'s tally each
    /// had their own copy with byte-identical bodies, and `open_record`'s justified
    /// withholding an `Ord` impl on the grounds that it "would export **this file's**
    /// sorting convention to every other consumer" — which the STR copy refutes: two
    /// independent producers wanted the identical convention, so it is the type's, not
    /// either file's (Milestone A review, Mi11). One copy means C4 rewrites this order
    /// once when the payload becomes a set, rather than twice with the two free to
    /// disagree while both still compile.
    ///
    /// The order it produces is the pre-reshape one — complete, then left-flush partials,
    /// then right-flush — without naming a side: a left-flush run has
    /// `offset_in_locus == 0` and so sorts ahead of every right-flush one, whose offset is
    /// `locus_len - covered`. Pinned by
    /// `witness_order_ranks_partials_by_offset_before_length`, which is the test that
    /// fails when the two components are exchanged.
    pub fn sort_key(self) -> (u8, u16, u16) {
        match self {
            Self::Complete => (0, 0, 0),
            Self::Partial {
                offset_in_locus,
                positions_covered,
            } => (1, offset_in_locus, positions_covered),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::collections::BTreeSet;
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};

    use proptest::prelude::*;

    fn hash_of(positions: &WitnessedLocusPositions) -> u64 {
        let mut hasher = DefaultHasher::new();
        positions.hash(&mut hasher);
        hasher.finish()
    }

    fn runs_of(positions: &WitnessedLocusPositions) -> Vec<(u16, u16)> {
        positions.runs().collect()
    }

    /// **Touching and overlapping runs are one run** — the merge half of the invariant.
    ///
    /// A gap of even one position is a real gap and survives: that is the distinction the
    /// whole type exists to carry, so it is asserted beside the merge rather than trusted.
    #[test]
    fn adjacent_and_overlapping_runs_merge_and_a_gap_survives() {
        let touching =
            WitnessedLocusPositions::from_half_open_runs([(0, 4), (4, 9)]).expect("non-empty");
        assert_eq!(
            runs_of(&touching),
            vec![(0, 9)],
            "runs that touch are one run"
        );

        let overlapping =
            WitnessedLocusPositions::from_half_open_runs([(0, 6), (3, 9)]).expect("non-empty");
        assert_eq!(runs_of(&overlapping), vec![(0, 9)]);

        let contained =
            WitnessedLocusPositions::from_half_open_runs([(0, 9), (3, 5)]).expect("non-empty");
        assert_eq!(
            runs_of(&contained),
            vec![(0, 9)],
            "a run inside another adds no positions and must not survive as a second run",
        );

        let gapped =
            WitnessedLocusPositions::from_half_open_runs([(0, 4), (5, 9)]).expect("non-empty");
        assert_eq!(
            runs_of(&gapped),
            vec![(0, 4), (5, 9)],
            "one unwitnessed position at 4 is a hole, and a hole is the point",
        );
    }

    /// **`LocusLen` saturates, and until this test nothing said so** — replacing the
    /// saturation with `positions as u16` left all 287 tests green (Milestone B review).
    ///
    /// The truncation is not a rounding error: a 65,536-position region would become
    /// `LocusLen(0)`, on which *every* witness is flush-right and zero-long, so the depth
    /// derivation reads zero everywhere and nothing panics. Saturating instead keeps the
    /// derivation total, which is the documented choice — a locus that long is a satellite
    /// the caller does not handle.
    #[test]
    fn locus_len_saturates_rather_than_truncating() {
        assert_eq!(LocusLen::from_positions(9).get(), 9);
        assert_eq!(
            LocusLen::from_positions(u64::from(u16::MAX)).get(),
            u16::MAX
        );
        assert_eq!(
            LocusLen::from_positions(u64::from(u16::MAX) + 1).get(),
            u16::MAX,
            "one past the top saturates; a truncating cast would give 0",
        );
        assert_eq!(
            LocusLen::from_positions(u64::MAX).get(),
            u16::MAX,
            "and so does anything above it",
        );
    }

    /// `of_region` is `from_positions` over the region's own length — the one source a
    /// consumer should use, so it is asserted to agree rather than assumed to.
    #[test]
    fn of_region_is_the_regions_length() {
        let region = GenomeRegion {
            contig: crate::ng::types::ContigId(0),
            start: crate::ng::types::Position(10),
            end: crate::ng::types::Position(20),
        };
        assert_eq!(region.len(), 11, "1-based, inclusive of both ends");
        assert_eq!(LocusLen::of_region(region), LocusLen::from_positions(11));
    }

    /// **A set nothing witnessed is not a set** — the "never empty" half of the invariant,
    /// and the run of no positions both `RefSpan` and `ReadWitness` document as not
    /// existing.
    #[test]
    fn an_empty_input_or_an_empty_run_is_rejected_rather_than_dropped() {
        assert!(WitnessedLocusPositions::from_half_open_runs([]).is_none());
        assert!(
            WitnessedLocusPositions::from_half_open_runs([(4, 4)]).is_none(),
            "a half-open run with start == end covers nothing",
        );
        assert!(
            WitnessedLocusPositions::from_half_open_runs([(0, 4), (7, 7)]).is_none(),
            "an empty run among good ones is rejected, not silently dropped — the caller \
             built a set it did not mean",
        );
        assert!(
            WitnessedLocusPositions::from_half_open_runs([(9, 4)]).is_none(),
            "and a reversed run is rejected too — `start >= end` covers both, which only a \
             test with `start > end` actually says",
        );
        assert!(WitnessedLocusPositions::one_run_from_offset_and_length(3, 0).is_none());
        assert!(
            WitnessedLocusPositions::one_run_from_offset_and_length(u16::MAX, 1).is_none(),
            "a run ending past u16::MAX fails rather than clamping (spec §3.4)",
        );
        assert!(
            WitnessedLocusPositions::one_run_from_offset_and_length(60_000, 10_000).is_none(),
            "and an end that overflows mid-range fails too. **This is the case that \
             discriminates**: at `u16::MAX + 1` the wrapped, saturated and checked sums all \
             land `<= start`, so `start >= end` rejects the run whatever the addition does, \
             and replacing `checked_add` with `saturating_add` left all 287 tests green \
             (Milestone B review). Here saturation would return a set of 5,535 positions — \
             the silent shortening the doc forbids",
        );
    }

    /// `one_run_from_offset_and_length` and `from_half_open_runs` describe the same set
    /// through the two spellings a caller must choose between. A reader meeting both in one
    /// type is entitled to know they agree.
    #[test]
    fn one_run_is_the_half_open_pair_it_looks_like() {
        assert_eq!(
            WitnessedLocusPositions::one_run_from_offset_and_length(4, 5),
            WitnessedLocusPositions::from_half_open_runs([(4, 9)]),
        );
    }

    /// **Positions, not span.** The sum skips the hole, which is the number
    /// `num_obs_along_locus` needs and the one a first-start-to-last-end subtraction gets
    /// wrong.
    #[test]
    fn positions_covered_counts_the_runs_and_not_the_gap_between_them() {
        let gapped =
            WitnessedLocusPositions::from_half_open_runs([(0, 3), (7, 9)]).expect("non-empty");
        assert_eq!(gapped.positions_covered(), 5);
        assert_eq!(
            WitnessedLocusPositions::from_half_open_runs([(0, 9)])
                .expect("non-empty")
                .positions_covered(),
            9,
            "and the same span with nothing missing counts every position",
        );
    }

    /// The flush predicates read off the first and last run, and mean what they meant on
    /// the one-run form: a prefix constraint, and a suffix constraint.
    #[test]
    fn flushness_is_read_off_the_first_and_last_run() {
        let len = LocusLen::from_positions(9);

        let whole = WitnessedLocusPositions::from_half_open_runs([(0, 9)]).expect("non-empty");
        assert!(whole.is_flush_left() && whole.is_flush_right(len));

        let interior =
            WitnessedLocusPositions::from_half_open_runs([(2, 3), (5, 7)]).expect("non-empty");
        assert!(!interior.is_flush_left() && !interior.is_flush_right(len));

        let both_ends_holed =
            WitnessedLocusPositions::from_half_open_runs([(0, 2), (7, 9)]).expect("non-empty");
        assert!(
            both_ends_holed.is_flush_left() && both_ends_holed.is_flush_right(len),
            "a read blind only in the middle is constrained from both sides",
        );
    }

    /// **The one-run predicate's addition must not wrap**, which nothing said until the
    /// Milestone B review swapped `saturating_add` for `wrapping_add` and watched all 287
    /// tests pass. `Partial`'s fields are public and its clamping is a convention rather
    /// than a type invariant (see the variant's own doc), so an over-long run is reachable
    /// without going through a constructor — and under wrapping it reports *not* flush
    /// right, which reads as "this read is missing the locus's tail" when it covered
    /// everything.
    #[test]
    fn a_run_whose_end_would_overflow_is_still_flush_right() {
        let over_long = ReadWitness::Partial {
            offset_in_locus: u16::MAX - 1,
            positions_covered: 8,
        };
        assert!(
            over_long.is_flush_right(LocusLen::from_positions(u64::from(u16::MAX))),
            "the reach saturates at u16::MAX and still reaches the border",
        );
    }

    /// **The set is one representation whether it fits inline or spills to the heap.**
    ///
    /// The encoding keeps two runs inline; a third spills. If equality could see the
    /// difference, a witness that arrived as three runs and merged down to two would stop
    /// matching the same two runs built directly — the exact silent split this type exists
    /// to prevent, on the exact path (a merge) that produces it.
    #[test]
    fn a_set_that_merged_down_to_two_runs_equals_the_same_two_built_directly() {
        let merged =
            WitnessedLocusPositions::from_half_open_runs([(0, 2), (5, 7), (1, 3), (6, 9), (2, 4)])
                .expect("non-empty");
        let direct =
            WitnessedLocusPositions::from_half_open_runs([(0, 4), (5, 9)]).expect("non-empty");
        assert_eq!(merged, direct);
        assert_eq!(hash_of(&merged), hash_of(&direct));
    }

    proptest! {
        /// **Two sets built from different orderings of the same runs compare and hash
        /// equal** — the property the whole design rests on, over generated multisets
        /// rather than a fixture.
        ///
        /// This is the guard the plan names for B2, and it is a *property* test because the
        /// failure is silent: an observation's identity is `(bases, read_witness,
        /// read_group)`, so a second representation of one set does not error, it splits one
        /// observation into two and inflates the count (spec §3.3).
        ///
        /// **Mutation-checked in three directions**: delete `sort_unstable` and this fails
        /// on a shuffled input; stop merging runs that touch and it fails on those; and
        /// truncate a containing run to the contained one — `open_end.max(end)` → `end` —
        /// and the *position-set* assertion below fails. A fixture would have caught
        /// neither of the first two reliably, because both need an input order the fixture
        /// happens not to use.
        ///
        /// **The position-set assertion is not decoration.** The first version of this test
        /// asserted only the canonical *shape* — sorted, disjoint, non-empty, and equal
        /// across orderings — and the Milestone B review showed that a normaliser which
        /// *loses witnessed positions* satisfies every one of them: sorting happens before
        /// merging, so both orderings lose the same positions and order-independence still
        /// holds. Shape is what the representation looks like; the positions are what it
        /// means.
        #[test]
        fn a_set_is_the_same_set_whatever_order_its_runs_arrive_in(
            runs in prop::collection::vec((0u16..40, 1u16..8), 1..6),
            permutation in prop::collection::vec(any::<u8>(), 6),
        ) {
            let half_open: Vec<(u16, u16)> = runs
                .iter()
                .map(|(start, len)| (*start, start + len))
                .collect();
            // A cheap deterministic shuffle: sort the same runs by a generated key, so the
            // second build sees an order the first did not.
            let mut shuffled: Vec<(u16, u16)> = half_open.clone();
            shuffled.sort_by_key(|run| {
                permutation[(usize::from(run.0) + usize::from(run.1)) % permutation.len()]
            });

            let first = WitnessedLocusPositions::from_half_open_runs(half_open.clone())
                .expect("at least one run");
            let second = WitnessedLocusPositions::from_half_open_runs(shuffled)
                .expect("at least one run");

            prop_assert_eq!(&first, &second, "one set, two orderings, one representation");
            prop_assert_eq!(hash_of(&first), hash_of(&second));

            // The canonical form itself, stated rather than assumed: sorted, non-empty,
            // and separated by at least one unwitnessed position.
            let canonical: Vec<(u16, u16)> = first.runs().collect();
            for pair in canonical.windows(2) {
                prop_assert!(
                    pair[0].1 < pair[1].0,
                    "runs {:?} touch or overlap, so they are one run",
                    pair,
                );
            }
            prop_assert!(canonical.iter().all(|(start, end)| start < end));

            // **And it is the same *positions*, not merely the same shape.** Every position
            // the input claimed, and no others.
            let claimed: BTreeSet<u16> = half_open
                .iter()
                .flat_map(|(start, end)| *start..*end)
                .collect();
            let kept: BTreeSet<u16> = canonical
                .iter()
                .flat_map(|(start, end)| *start..*end)
                .collect();
            prop_assert_eq!(
                kept,
                claimed,
                "canonicalising changed which positions are witnessed",
            );
            prop_assert_eq!(
                first.positions_covered(),
                u32::try_from(canonical.iter().map(|(s, e)| usize::from(*e - *s)).sum::<usize>())
                    .expect("a locus is at most u16::MAX positions"),
            );
        }
    }

    /// **The sort key's two components, told apart** — and until this test nothing could
    /// tell them apart.
    ///
    /// [`sort_key`](ReadWitness::sort_key) is what makes the emitted order a function of
    /// an observation's own identity rather than of read arrival order. The Milestone A
    /// review mutated it — exchanging `offset_in_locus` and `positions_covered` in the
    /// returned key — and **all 275 tests stayed green**, while a `panic!` in the same arm
    /// failed four of them: the arm runs, and nothing asserted what it returned. The test
    /// named for the job, `parity::the_projection_orders_observations_as_the_walk_does`,
    /// sorts both sides with *this same key*, so it cannot fail on any change to it.
    ///
    /// The discriminating input is two runs whose components vary in **opposite**
    /// directions, which no fixture in the suite produced: `{0, 9}` against `{4, 2}`. If
    /// length outranked offset, the long left-flush run would sort second.
    #[test]
    fn witness_order_ranks_partials_by_offset_before_length() {
        assert!(
            ReadWitness::Complete.sort_key()
                < ReadWitness::Partial {
                    offset_in_locus: 0,
                    positions_covered: 1,
                }
                .sort_key(),
            "a complete witness sorts ahead of every run",
        );
        assert!(
            ReadWitness::Partial {
                offset_in_locus: 0,
                positions_covered: 9,
            }
            .sort_key()
                < ReadWitness::Partial {
                    offset_in_locus: 4,
                    positions_covered: 2,
                }
                .sort_key(),
            "where a run starts outranks how far it runs, so a left-flush run precedes a \
             right-flush one whatever their lengths",
        );
    }

    /// **The constructors place the run against their own border**, and the offset is
    /// derived from the *clamped* reach, never the raw one.
    ///
    /// The right case is the one that matters: `from_right(4, 10)` must start at 6. Deriving
    /// the offset before clamping would wrap on an over-long reach, and a saturating
    /// subtraction would quietly relabel a right-anchored read as left-anchored — the
    /// silent-depth failure this step was flagged for.
    #[test]
    fn from_right_places_the_run_against_the_right_border() {
        assert_eq!(
            ReadWitness::from_left(4, LocusLen::from_positions(10)),
            ReadWitness::Partial {
                offset_in_locus: 0,
                positions_covered: 4
            }
        );
        assert_eq!(
            ReadWitness::from_right(4, LocusLen::from_positions(10)),
            ReadWitness::Partial {
                offset_in_locus: 6,
                positions_covered: 4
            }
        );
    }

    /// **Once the reach covers the whole locus the two constructors agree** — a read that
    /// witnessed every position is constrained from neither border, so there is one run and
    /// not two.
    ///
    /// Stated as a test because it is a real behaviour change, not a curiosity: on the STR
    /// path an expanded allele can give a reach longer than the reference tract, and a
    /// left-anchored and a right-anchored read of the *same bases* then land in the **same
    /// tally bucket** and merge into one observation — where `PartialLeft(n)` and
    /// `PartialRight(n)` kept them apart. See `ssr::tally::tests::an_expanded_allele_merges_the_two_sides`.
    #[test]
    fn from_left_and_from_right_agree_once_the_reach_covers_the_whole_locus() {
        assert_eq!(
            ReadWitness::from_left(9, LocusLen::from_positions(3)),
            ReadWitness::Partial {
                offset_in_locus: 0,
                positions_covered: 3
            }
        );
        assert_eq!(
            ReadWitness::from_right(9, LocusLen::from_positions(3)),
            ReadWitness::from_left(9, LocusLen::from_positions(3))
        );
    }

    /// The flushness predicates — the **entire** surviving representation of
    /// prefix-versus-suffix, since the reshape dropped the side-tagged variants.
    ///
    /// Both `Complete` arms are asserted here because no call site reaches them: every
    /// `witness_label` matches `Complete` first, so inverting either arm to `false` would
    /// otherwise change nothing anywhere in the tree.
    #[test]
    fn flushness_is_derived_from_where_the_run_sits() {
        assert!(ReadWitness::Complete.is_flush_left());
        assert!(ReadWitness::Complete.is_flush_right(LocusLen::from_positions(10)));

        let left = ReadWitness::from_left(4, LocusLen::from_positions(10));
        assert!(left.is_flush_left(), "a prefix constraint");
        assert!(!left.is_flush_right(LocusLen::from_positions(10)));

        let right = ReadWitness::from_right(4, LocusLen::from_positions(10));
        assert!(!right.is_flush_left());
        assert!(
            right.is_flush_right(LocusLen::from_positions(10)),
            "a suffix constraint"
        );

        // An interior run — flush with neither border. The STR path cannot mint one, but the
        // predicates are shared and the generic path will.
        let interior = ReadWitness::Partial {
            offset_in_locus: 3,
            positions_covered: 4,
        };
        assert!(!interior.is_flush_left());
        assert!(!interior.is_flush_right(LocusLen::from_positions(10)));
    }
}
