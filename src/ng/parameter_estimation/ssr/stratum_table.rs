//! One stratum's evidence: every distinct locus shape and how many loci had it.
//!
//! **An entry is a locus, not a read**, and that is the design's central choice rather than
//! a storage detail. A read carries no genotype — it drew one of its locus's alleles and
//! then slipped — so a tally that pools reads across loci holds the allele spectrum
//! convolved with the slippage kernel, and recovering the kernel from that means undoing a
//! convolution with both halves unknown. Measured, the fitted slippage level then moves
//! **333-fold depending only on where the search starts**; keyed by locus the same fit is
//! exactly unbiased (`spec/parameter_prepass_ssr.md` §4.1).
//!
//! **Two words, two things.** A *locus shape* is how one locus's reads fell across the offset
//! buckets; an *entry* is one shape together with how many loci had it. This file holds the
//! first; the table of entries lands in Milestone B2.

use crate::ng::types::DomainError;

use super::{MAX_LOCUS_READS, OFFSET_BUCKETS, OffsetBucket, WholeRepeatOffset, bucket_of};

/// One locus's reads, laid out across the offset buckets: how many of them showed each
/// whole-repeat offset from the reference tract length, plus how many showed a length that is
/// not a whole number of copies.
///
/// **This is what the table is keyed on**, so two loci whose reads fell the same way become one
/// entry — one shape with a count of two — and which loci they were is never asked again.
///
/// The reads in the buckets plus the reads in the guard are the locus's depth as this unit saw
/// it, and nothing stores that sum: [`depth`](Self::depth) recomputes it. **The depth is exact
/// rather than binned, and that follows from the cap rather than from a separate decision** —
/// [`MAX_LOCUS_READS`] bounds it at a dozen values, so a ladder over it would save nothing. (If
/// the cap ever rises far, the generic path's measurement applies and is not free: its depth
/// ladder turned out to be a *correctness* parameter, sixteen bins costing 0.55 rungs of its
/// error-rate ladder — that ladder runs Phred 10 to 50 in 161 rungs, so a rung is a
/// quarter-Phred and 0.55 of one is 0.14 Phred — where twenty bins cost 0.05.)
///
/// **Ordered and hashable, which is the whole of the determinism requirement**: it keys a
/// `BTreeMap`, so the order two entries are visited in is a property of their contents and not
/// of when they arrived. Every fit is a floating-point sum over the entries, and
/// floating-point addition is not associative — a table walked in a different order between
/// two runs would let a fitted level wobble on identical data.
///
/// **A guard, not a parameter.** The non-whole-repeat count is modelled as an independent
/// per-read outcome, so the likelihood splits exactly into *how many reads showed a
/// non-whole-repeat length* times *how the rest fell across the offsets*: nothing about
/// slippage is estimated from it and nothing about it disturbs the slippage parameters
/// (`spec/parameter_prepass_ssr.md` §4.1). It is carried because §5's diagnostic is the share
/// of off-reference reads that land in it.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct LocusShape {
    /// Reads at each bucket, in bucket order; index `OFFSET_HALF_RANGE` is the reference
    /// length.
    reads_by_bucket: [u8; OFFSET_BUCKETS],
    /// Reads differing from the reference tract length by something other than a whole number
    /// of motif copies.
    reads_not_whole_repeat: u8,
}

/// **The cap and the counter width are one decision, and this is where it is held.** A shape
/// counts its reads in `u8`. Raise [`MAX_LOCUS_READS`] past 255 and [`LocusShape::try_new`]
/// admits a locus it cannot store: the narrowing there then fails and panics — loudly, and in
/// a release build that aborts rather than unwinds — on the first locus deep enough to reach
/// it. This assertion turns that into an `error[E0080]` at the edit that caused it.
const _: () = assert!(
    MAX_LOCUS_READS <= u8::MAX as u32,
    "a locus shape counts its reads in u8, so the read cap must fit in one"
);

impl LocusShape {
    /// The only constructor, and it **takes a wider integer than the shape stores**: the counts
    /// arrive as `u32`, the width [`MAX_LOCUS_READS`] itself is written in and the width a
    /// caller tallying a deep locus's reads will hold.
    ///
    /// A constructor as narrow as the storage would make `try_new(tally as u8)` the natural
    /// call, under which 260 reads arrive as 4 and validate as an ordinary shallow locus. Taking
    /// the wider type and rejecting instead means a caller that forgot to subsample gets a
    /// `Result`. (`DomainError::SsrPeriod` carries a `usize` against a `u8` period for the same
    /// *shape* of reason, though not the same one: there the width is what every producer of a
    /// period already holds.)
    ///
    /// Rejects a shape holding no reads at all as well as one above the cap. **An empty shape
    /// is not a harmless zero:** it would enter the table as a locus, so it would count
    /// towards `MIN_LOCI_TO_FIT` and towards every "how many loci stood behind this fit"
    /// number, while contributing a likelihood of exactly one to every candidate. A locus
    /// whose reads all failed to witness the tract is a locus with no evidence, and the
    /// caller's answer is to enter nothing rather than to enter this.
    pub fn try_new(
        reads_by_bucket: [u32; OFFSET_BUCKETS],
        reads_not_whole_repeat: u32,
    ) -> Result<Self, DomainError> {
        // Summed in `u64` because the arguments are unchecked: nine `u32`s and a tenth
        // overflow a `u32` long before they overflow this, so the total that reaches the
        // comparison is the total the caller offered rather than a wrapped remainder.
        let reads: u64 = reads_by_bucket
            .into_iter()
            .chain([reads_not_whole_repeat])
            .map(u64::from)
            .sum();

        if reads == 0 || reads > u64::from(MAX_LOCUS_READS) {
            return Err(DomainError::SsrLocusShapeReads {
                reads,
                cap: MAX_LOCUS_READS,
            });
        }

        // PANIC-FREE: every count is one term of the total, so none exceeds it; the check above
        // holds that total at or below `MAX_LOCUS_READS`, and the compile-time assertion above
        // holds `MAX_LOCUS_READS` inside a byte. Both narrowings rest on this one argument.
        let narrow = |count: u32| {
            u8::try_from(count).expect("no count can exceed the checked total of all counts")
        };
        Ok(Self {
            reads_by_bucket: reads_by_bucket.map(narrow),
            reads_not_whole_repeat: narrow(reads_not_whole_repeat),
        })
    }

    /// Reads at each bucket, in bucket order.
    ///
    /// Handed out as the whole array rather than one bucket at a time because that is how the
    /// fit reads it: a shape's contribution is a sum of `count · ln p` over the buckets, zipped
    /// against the per-bucket probabilities of one genotype. **This is the one place a count
    /// comes back in the width it is stored in** — the array is what the fit converts to `f64`
    /// element by element, so widening it here would build a second array for nothing.
    #[inline]
    #[must_use]
    pub fn reads_by_bucket(self) -> [u8; OFFSET_BUCKETS] {
        self.reads_by_bucket
    }

    /// Reads recorded in one bucket. Total: any [`OffsetBucket`] indexes a shape, which is what
    /// that type's private field buys.
    #[inline]
    #[must_use]
    pub fn reads_in(self, bucket: OffsetBucket) -> u32 {
        u32::from(self.reads_by_bucket[bucket.index()])
    }

    /// Reads that showed a length differing from the reference tract's by something other than
    /// a whole number of motif copies.
    #[inline]
    #[must_use]
    pub fn reads_not_whole_repeat(self) -> u32 {
        u32::from(self.reads_not_whole_repeat)
    }

    /// Reads whose length differs from the reference tract's by a whole number of motif copies,
    /// so they landed in an offset bucket — **the depth the fit scores over**, which is the
    /// locus's depth less the reads the guard holds. The design calls it the *scored depth*
    /// (`arch/parameter_prepass_ssr.md` §2.4).
    ///
    /// It can be zero at a shape that is not empty: a locus every one of whose reads showed a
    /// non-whole-repeat length is a real observation, and the length half of its likelihood is
    /// an empty product.
    #[inline]
    #[must_use]
    pub fn whole_repeat_depth(self) -> u32 {
        self.reads_by_bucket.iter().copied().map(u32::from).sum()
    }

    /// Every read this locus contributed, guard included. At least one, and at most
    /// [`MAX_LOCUS_READS`].
    #[inline]
    #[must_use]
    pub fn depth(self) -> u32 {
        self.whole_repeat_depth() + self.reads_not_whole_repeat()
    }

    /// Reads that showed a length other than the reference tract's — whole-repeat or not.
    ///
    /// **The denominator of §5's guard diagnostic**, and the numerator of nothing. Summed over
    /// every bucket but the reference one, plus the guard, so it is total by construction:
    /// written as a subtraction from the depth it would be a wrapping one if the reference
    /// bucket ever exceeded it, and a release build leaves `overflow-checks` off, so that would
    /// report a share near four billion rather than failing.
    #[inline]
    #[must_use]
    pub fn reads_off_reference(self) -> u32 {
        let reference = bucket_of(WholeRepeatOffset(0)).index();
        let off_reference: u32 = self
            .reads_by_bucket
            .iter()
            .enumerate()
            .filter(|&(bucket, _)| bucket != reference)
            .map(|(_, &reads)| u32::from(reads))
            .sum();
        off_reference + self.reads_not_whole_repeat()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::ssr::OFFSET_HALF_RANGE;

    /// A shape's reads, written the way a test reads them: the offsets that hold reads, and how
    /// many. Everything else is zero. Two offsets that saturate into the same end bucket add,
    /// which is what the bucket does.
    fn shape_with_reads(
        offsets: &[(i8, u32)],
        reads_not_whole_repeat: u32,
    ) -> Result<LocusShape, DomainError> {
        let mut reads_by_bucket = [0; OFFSET_BUCKETS];
        for &(offset, reads) in offsets {
            reads_by_bucket[bucket_of(WholeRepeatOffset(offset)).index()] += reads;
        }
        LocusShape::try_new(reads_by_bucket, reads_not_whole_repeat)
    }

    /// The invariant the whole type rests on: nothing stores the depth, so the reads a shape
    /// reports are exactly the reads it was built from — split into the ones the fit scores and
    /// the ones the guard holds.
    #[test]
    fn a_shape_reports_the_reads_it_was_built_from_split_at_the_guard() {
        let shape =
            shape_with_reads(&[(0, 4), (-1, 2), (2, 1)], 3).expect("ten reads, inside the cap");

        assert_eq!(shape.whole_repeat_depth(), 7);
        assert_eq!(shape.reads_not_whole_repeat(), 3);
        assert_eq!(shape.depth(), 10);
        // Against the literal the fixture was built from, so the two sides are independent: an
        // assertion summing `reads_by_bucket()` and comparing it with `whole_repeat_depth()`
        // compares that function with a copy of itself and cannot fail.
        assert_eq!(
            shape.reads_by_bucket(),
            [0, 0, 0, 2, 4, 0, 1, 0, 0],
            "each read sits in the bucket of the offset it showed"
        );
    }

    /// Every read is recorded where it was put — a shape that transposed two buckets would
    /// report a stratum losing repeats as one gaining them, which is the parameter this whole
    /// path exists to protect.
    #[test]
    fn each_read_is_recorded_in_the_bucket_of_the_offset_it_showed() {
        let shape = shape_with_reads(&[(0, 5), (-2, 3), (1, 1)], 0).expect("nine reads");

        assert_eq!(shape.reads_in(bucket_of(WholeRepeatOffset(0))), 5);
        assert_eq!(shape.reads_in(bucket_of(WholeRepeatOffset(-2))), 3);
        assert_eq!(shape.reads_in(bucket_of(WholeRepeatOffset(1))), 1);
        assert_eq!(shape.reads_in(bucket_of(WholeRepeatOffset(3))), 0);
    }

    /// **Every bucket holds its own reads, and the two accessors agree on all nine.** The fit
    /// zips [`LocusShape::reads_by_bucket`] against a per-bucket probability while the
    /// diagnostics reach for one bucket through [`LocusShape::reads_in`], so a bucket the two
    /// read differently is a fit scoring one thing and a report naming another. Four buckets
    /// are named by the fixtures above; this reaches the other five.
    #[test]
    fn every_bucket_round_trips_through_both_accessors() {
        let mut reads_by_bucket = [0u32; OFFSET_BUCKETS];
        for (bucket, reads) in reads_by_bucket.iter_mut().enumerate() {
            *reads = u32::from(bucket % 2 == 0);
        }
        let shape =
            LocusShape::try_new(reads_by_bucket, 2).expect("five scored reads and two guard");

        for offset in -OFFSET_HALF_RANGE..=OFFSET_HALF_RANGE {
            let bucket = bucket_of(WholeRepeatOffset(offset));
            assert_eq!(
                shape.reads_in(bucket),
                reads_by_bucket[bucket.index()],
                "offset {offset} did not come back from the bucket it went into"
            );
            assert_eq!(
                u32::from(shape.reads_by_bucket()[bucket.index()]),
                shape.reads_in(bucket)
            );
        }
        assert_eq!(shape.whole_repeat_depth(), 5);
        assert_eq!(shape.depth(), 7);
    }

    /// **Reads that saturated into an end bucket are ordinary reads**, at both ends: stored,
    /// read back through both accessors, and counted in both depths.
    ///
    /// The end buckets are where a locus carrying a long or a short allele puts its reads —
    /// 0.89% of reads across GRCh38's typed tracts — so a bucket dropped here is a locus whose
    /// evidence quietly shrinks rather than one that fails. Every other fixture in this file
    /// sits at an offset the range holds exactly, which is why this one exists.
    #[test]
    fn reads_that_saturated_into_an_end_bucket_are_stored_read_back_and_counted() {
        let far_long = shape_with_reads(&[(0, 2), (7, 3)], 1).expect("six reads");

        assert_eq!(far_long.reads_in(bucket_of(WholeRepeatOffset(7))), 3);
        assert_eq!(far_long.reads_by_bucket()[OFFSET_BUCKETS - 1], 3);
        assert_eq!(far_long.whole_repeat_depth(), 5);
        assert_eq!(far_long.depth(), 6);
        assert_eq!(far_long.reads_off_reference(), 4);

        let far_short = shape_with_reads(&[(0, 2), (-9, 3)], 1).expect("six reads");

        assert_eq!(far_short.reads_in(bucket_of(WholeRepeatOffset(-9))), 3);
        assert_eq!(far_short.reads_by_bucket()[0], 3);
        assert_eq!(far_short.whole_repeat_depth(), 5);
        assert_eq!(far_short.depth(), 6);
        assert_eq!(far_short.reads_off_reference(), 4);
    }

    /// **The cap is on the locus's reads and not on the ones the fit scores**, so the guard
    /// counts towards it: nine scored reads plus four guard reads is a locus of thirteen, above
    /// a cap of twelve, whichever channel they land in.
    #[test]
    fn a_shape_above_the_read_cap_cannot_be_built_and_the_guard_counts_towards_it() {
        let at_the_cap =
            shape_with_reads(&[(0, 8), (-1, 4)], 0).expect("twelve reads is the cap, not past it");
        assert_eq!(
            at_the_cap.depth(),
            MAX_LOCUS_READS,
            "a shape at the cap reports the cap"
        );

        for (offsets, guard, reads) in [
            (&[(0i8, 13u32)][..], 0u32, 13u64),
            (&[(0, 9)][..], 4, 13),
            (&[(0, 6), (-1, 7)][..], 0, 13),
        ] {
            assert_eq!(
                shape_with_reads(offsets, guard),
                Err(DomainError::SsrLocusShapeReads {
                    reads,
                    cap: MAX_LOCUS_READS
                }),
                "{reads} reads, of which {guard} in the guard"
            );
        }
    }

    /// **One read is a locus with evidence**, and it is the lower edge of the accepted range.
    /// A floor set anywhere above one would drop the shallowest loci silently: the caller's
    /// answer to a rejected shape is to enter nothing, so those loci would simply not appear in
    /// the stratum's count or in what `MIN_LOCI_TO_FIT` is measured against.
    #[test]
    fn a_shape_of_a_single_read_is_accepted() {
        let one_scored = shape_with_reads(&[(0, 1)], 0).expect("one read is evidence");
        assert_eq!(one_scored.depth(), 1);
        assert_eq!(one_scored.whole_repeat_depth(), 1);

        let one_guard = shape_with_reads(&[], 1).expect("one guard read is evidence");
        assert_eq!(one_guard.depth(), 1);
        assert_eq!(one_guard.whole_repeat_depth(), 0);

        let two = shape_with_reads(&[(0, 1), (1, 1)], 0).expect("two reads");
        assert_eq!(two.depth(), 2);
    }

    /// **A count above a byte is rejected rather than truncated**, which is the reason the
    /// constructor takes a wider integer than it stores. Narrowed first, 260 reads at the
    /// reference length would arrive as 4 and validate as an ordinary shallow locus — a depth
    /// wrong 65-fold with nothing to notice it.
    ///
    /// What this pins is the **signature**: it would not compile against a `try_new([u8; N], u8)`.
    /// The rejection itself comes from the ordinary cap comparison, the same one that rejects
    /// thirteen — with the cap inside a byte there is no separate width test to reach.
    #[test]
    fn a_count_wider_than_the_storage_is_rejected_rather_than_wrapped() {
        assert_eq!(
            shape_with_reads(&[(0, 260)], 0),
            Err(DomainError::SsrLocusShapeReads {
                reads: 260,
                cap: MAX_LOCUS_READS
            })
        );
        assert_eq!(
            shape_with_reads(&[], 300),
            Err(DomainError::SsrLocusShapeReads {
                reads: 300,
                cap: MAX_LOCUS_READS
            })
        );
    }

    /// The offered total is summed wide enough to be reported, so a caller handed a message
    /// sees the number it offered rather than a wrapped remainder of it. The fixture offers ten
    /// values of `u32::MAX` — nine buckets and the guard — which is ten times what a `u32` sum
    /// could hold.
    #[test]
    fn the_rejection_names_the_total_the_caller_offered_however_large() {
        let every_bucket_full = LocusShape::try_new([u32::MAX; OFFSET_BUCKETS], u32::MAX);

        assert_eq!(
            every_bucket_full,
            Err(DomainError::SsrLocusShapeReads {
                reads: 10 * u64::from(u32::MAX),
                cap: MAX_LOCUS_READS
            })
        );
    }

    /// **An empty shape is refused**, because it would enter the table as a locus: it would
    /// count towards `MIN_LOCI_TO_FIT` and towards every "loci behind this fit" number while
    /// scoring exactly one under every candidate. A locus none of whose reads witnessed the
    /// tract is a locus to leave out, not one to enter empty.
    #[test]
    fn a_shape_holding_no_reads_at_all_is_refused() {
        assert_eq!(
            shape_with_reads(&[], 0),
            Err(DomainError::SsrLocusShapeReads {
                reads: 0,
                cap: MAX_LOCUS_READS
            })
        );
    }

    /// A locus every one of whose reads showed a non-whole-repeat length is a real observation
    /// and not an empty shape — the fit scores none of it, and the guard diagnostic is exactly
    /// what it is there to report.
    #[test]
    fn a_shape_whose_reads_are_all_non_whole_repeat_is_still_a_real_shape() {
        let all_guard = shape_with_reads(&[], 5).expect("five reads, all of them non-whole-repeat");

        assert_eq!(all_guard.whole_repeat_depth(), 0);
        assert_eq!(all_guard.depth(), 5);
        assert_eq!(all_guard.reads_off_reference(), 5);
    }

    /// The guard's reads differ from the reference length too, so they belong in the
    /// denominator of §5's share along with the whole-repeat ones. A shape counting only the
    /// off-reference *buckets* would put the guard in the numerator and not the denominator,
    /// and report shares above one.
    #[test]
    fn reads_off_reference_counts_the_guard_as_well_as_the_other_buckets() {
        let shape = shape_with_reads(&[(0, 6), (-1, 2)], 1).expect("nine reads");

        assert_eq!(shape.reads_off_reference(), 3);
        assert_eq!(
            shape_with_reads(&[(0, 9)], 0)
                .expect("nine reads at the reference")
                .reads_off_reference(),
            0
        );
    }

    /// **Two loci that looked alike are one key**, which is what makes the table collapse 1.73
    /// million tomato loci into a small object — and shapes that differ anywhere, the guard
    /// included, are two. A type whose equality ignored the guard would merge a clean locus
    /// with an interrupted one.
    #[test]
    fn shapes_are_equal_and_hash_alike_exactly_when_every_count_agrees() {
        use std::collections::HashSet;

        let clean = shape_with_reads(&[(0, 4)], 0).expect("four reads");
        let same = shape_with_reads(&[(0, 4)], 0).expect("four reads");
        let interrupted = shape_with_reads(&[(0, 4)], 1).expect("four reads and a guard read");
        let shifted = shape_with_reads(&[(-1, 4)], 0).expect("four reads, a repeat short");

        assert_eq!(clean, same);
        assert_ne!(clean, interrupted);
        assert_ne!(clean, shifted);

        let distinct: HashSet<LocusShape> =
            [clean, same, interrupted, shifted].into_iter().collect();
        assert_eq!(distinct.len(), 3, "equal shapes must hash alike");
    }

    /// **A key that ignored one bucket would merge two loci that observed different things**,
    /// and the merge is silent: the entry's count rises by one and the shape it is filed under
    /// is not the shape either locus had. The test above reaches two of the ten places a shape
    /// can differ; this sweeps all ten, through both the hashing the type derives and the
    /// `BTreeMap` ordering the table will key on.
    #[test]
    fn shapes_differing_in_any_single_bucket_are_distinct_keys() {
        use std::collections::{BTreeMap, HashSet};

        let mut shapes = Vec::new();
        for bucket in 0..OFFSET_BUCKETS {
            let mut reads_by_bucket = [1u32; OFFSET_BUCKETS];
            reads_by_bucket[bucket] = 2;
            shapes.push(LocusShape::try_new(reads_by_bucket, 0).expect("ten reads"));
        }
        shapes.push(LocusShape::try_new([1; OFFSET_BUCKETS], 0).expect("nine reads"));
        shapes.push(LocusShape::try_new([1; OFFSET_BUCKETS], 1).expect("ten reads"));

        let offered = shapes.len();
        assert_eq!(
            shapes.iter().copied().collect::<HashSet<_>>().len(),
            offered,
            "a bucket the key ignores would merge two loci that observed different things"
        );
        assert_eq!(
            shapes
                .iter()
                .map(|shape| (*shape, ()))
                .collect::<BTreeMap<_, _>>()
                .len(),
            offered,
            "a bucket the ordering ignores would merge two loci in the table"
        );
    }

    /// **The order separates every pair of distinct shapes**, which is what a `BTreeMap` key
    /// has to do. That the order is content-determined rather than arrival-determined is
    /// discharged by `Ord` being derived — nothing a test can hand `Ord` carries an arrival
    /// time — and what this fixture can fail on is an order that ties two distinct shapes:
    /// `sort` is stable, so a tie reproduces the input order and the two sorts differ. An `Ord`
    /// keyed on the depth alone, which is a perfectly good total order, ties three of these
    /// four.
    #[test]
    fn shapes_sort_the_same_way_whatever_order_they_arrived_in() {
        let shapes = [
            shape_with_reads(&[(0, 3)], 0).expect("three reads"),
            shape_with_reads(&[(0, 2), (1, 1)], 0).expect("three reads"),
            shape_with_reads(&[(-4, 3)], 0).expect("three reads"),
            shape_with_reads(&[(0, 3)], 1).expect("four reads"),
        ];

        let mut forwards = shapes;
        forwards.sort();
        let mut backwards = shapes;
        backwards.reverse();
        backwards.sort();

        assert_eq!(forwards, backwards);
        assert_eq!(
            *forwards.last().expect("four shapes"),
            shape_with_reads(&[(-4, 3)], 0).expect("three reads"),
            "comparison starts at bucket 0, so the only shape with reads there sorts last"
        );
    }
}
