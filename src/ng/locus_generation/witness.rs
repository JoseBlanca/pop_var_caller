//! What one read saw of one locus — the witness vocabulary.
//!
//! [`ReadWitness`] is the answer to "how much of this locus did this read witness", and
//! [`LocusLen`] is the axis it is measured on. They live in their own file rather than in
//! `mod.rs` because the invariant that matters here — one set of witnessed positions has
//! exactly one representation, so two reads that witnessed the same thing share one
//! observation — fails *silently* when it fails, and a private field needs one module
//! boundary to be private to (arch *Module home*).
//!
//! Re-exported from [`locus_generation`](super), so no consumer's import path names this
//! module.

use crate::ng::types::GenomeRegion;

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
