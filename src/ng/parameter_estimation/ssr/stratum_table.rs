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
//! buckets; an *entry* is one shape together with how many loci had it. This file holds both,
//! and the second is [`StratumTable`].
//!
//! Nothing here reads a locus: a shape and its two base counts arrive already made, from
//! Milestone C. What this file owns is that they store, merge and report exactly.

use std::collections::BTreeMap;

use crate::ng::types::{DomainError, ErrorRate};

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

/// How many bases of one locus's reads were compared against the tract they sit in, and how
/// many of those bases differed from it.
///
/// **The second thing a locus contributes**, beside its shape, and the whole of what the
/// per-base substitution rate is fitted from. The two counts are a sufficient statistic for it:
/// a read's mismatch count is binomial at that rate whatever length the read showed, so the
/// channel that scores *which offset a read landed at* and the channel that scores *how its
/// bases read* factorise exactly, and the rate is a division rather than an axis of any search
/// (`spec/parameter_prepass_ssr.md` §4.1).
///
/// **One type, where the architecture sketches two arguments.** Two same-typed counts side by
/// side are a transposition waiting to happen, and one is a subset of the other, so a swapped
/// pair is not obviously wrong at the call site. What a swapped pair does downstream depends on
/// what it is pooled with, and neither outcome is a report a reader could act on: alone it makes
/// the stratum's rate exceed one, which [`ErrorRate`] refuses, so the run dies inside
/// [`StratumTable::substitution_rate`]; pooled with well-formed loci it survives as a plausible
/// wrong rate that nothing downstream can question. Here the pair cannot be built the wrong way
/// round — mismatched bases outnumbering compared ones are refused — which is also what makes
/// the fitted rate a probability by construction rather than by hope.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct BaseComparison {
    bases_compared: u32,
    bases_mismatched: u32,
}

impl BaseComparison {
    /// The only constructor. More mismatched bases than compared ones is not a noisy locus: it
    /// is a caller that has swapped its two counts, or counted one of them wrongly.
    pub fn try_new(bases_compared: u32, bases_mismatched: u32) -> Result<Self, DomainError> {
        if bases_mismatched > bases_compared {
            return Err(DomainError::SsrBaseComparison {
                bases_compared,
                bases_mismatched,
            });
        }
        Ok(Self {
            bases_compared,
            bases_mismatched,
        })
    }

    /// Bases of this locus's reads that were compared against the tract.
    #[inline]
    #[must_use]
    pub fn bases_compared(self) -> u32 {
        self.bases_compared
    }

    /// Of those, the ones that differed. Never more than [`Self::bases_compared`].
    #[inline]
    #[must_use]
    pub fn bases_mismatched(self) -> u32 {
        self.bases_mismatched
    }
}

/// One entry of a stratum's table, as a fit reads it: a locus shape, and how many loci had it.
///
/// **Named fields rather than the `(LocusShape, u64)` the architecture sketches**, for the reason
/// the sibling path's `Cell` replaced its own sketched tuple: nothing in a tuple says that the
/// second member is a count of *loci* rather than of reads or of entries, and the fit re-walks
/// this list once per candidate in code that would otherwise read `.0` and `.1`. Nothing outside
/// this file consumes it yet, so this is the cheapest moment it will ever have.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct StratumEntry {
    /// How one locus's reads fell across the offset buckets.
    pub shape: LocusShape,
    /// How many loci in this stratum showed that shape.
    pub loci: u64,
}

/// One stratum's evidence, for one read group: every distinct locus shape and how many loci had
/// it, plus the two base counts the substitution rate is fitted from.
///
/// **Sparse rather than dense, and that is forced rather than chosen.** The sibling path's cell
/// space is a depth bin against an alternative-read count — 583 cells, worth holding densely.
/// Here an entry is a whole locus's split across ten buckets, so the space of possible shapes at
/// one depth is 220 at three reads a locus and 293,930 at twelve
/// (`arch/parameter_prepass_ssr.md` §2.2), of which only a small, data-dependent corner is ever
/// occupied: a whole tomato genome's 1.73 million STR loci made 70,305 entries — one entry per 25
/// loci — and HG002 at 300× made 12,727 for 29,811 loci, one per 2.3
/// (`spec/parameter_prepass_ssr.md` §4.1). Most loci at a clean tract are "every read at the
/// reference length", and what separates two entries is then mostly their depth.
///
/// **A `BTreeMap` and not a `HashMap`**, for the reason [`LocusShape`] is ordered: every fit is a
/// floating-point sum over the entries and floating-point addition is not associative, so the
/// walk order has to be a property of the contents rather than of a hash seed or of arrival.
///
/// **A locus covered by two read groups makes one entry in each of two tables, and that is
/// sound.** Each table's own distribution is correctly specified — the genotype is drawn once for
/// the locus and enters both through the same mixture — so the product over them is a composite
/// likelihood: the estimator stays consistent and what the split throws away is dependence
/// between the two, which is precision. The sibling path measured exactly that split and found it
/// unbiased in all 25 simulated samples tried (`spec/parameter_prepass_ssr.md` §4.1). What must
/// *not* be split is a locus's reads within one read group, which is what keying an entry by
/// locus exists to prevent.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct StratumTable {
    /// How many loci showed each shape. Never holds a zero count: an entry exists because a
    /// locus made it.
    loci_by_shape: BTreeMap<LocusShape, u32>,
    /// Every locus's [`BaseComparison`], added up.
    bases: PooledBases,
}

/// Every locus's [`BaseComparison`] in one stratum, added up.
///
/// **The same two counts as [`BaseComparison`], kept as one type for the same reason**: they are
/// the same width and one is a subset of the other, so a pair carried positionally through a
/// helper transposes without a compiler complaint.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
struct PooledBases {
    compared: u64,
    mismatched: u64,
}

impl PooledBases {
    /// Add one more locus's comparison, or another shard's pool.
    ///
    /// # Panics
    ///
    /// If either counter passes `u64`, which no genome reaches. Kept because the release profile
    /// leaves `overflow-checks` off, so an unchecked `+=` would wrap a stratum's rate to a number
    /// computed from a remainder — which nothing downstream could question.
    fn add(&mut self, more: Self) {
        self.compared = self.compared.checked_add(more.compared).unwrap_or_else(|| {
            panic!("bases compared in one stratum passed {}", u64::MAX);
        });
        self.mismatched = self
            .mismatched
            .checked_add(more.mismatched)
            .unwrap_or_else(|| {
                panic!("bases mismatched in one stratum passed {}", u64::MAX);
            });
    }
}

impl From<BaseComparison> for PooledBases {
    fn from(bases: BaseComparison) -> Self {
        Self {
            compared: u64::from(bases.bases_compared()),
            mismatched: u64::from(bases.bases_mismatched()),
        }
    }
}

impl StratumTable {
    /// Add one locus: its shape, and the bases its reads compared against the tract.
    ///
    /// Never allocates after the first locus of a given shape, which is what the measured
    /// collapse buys — at 0.43 entries a locus on HG002, 57 loci in 100 find an entry already
    /// there.
    ///
    /// # Panics
    ///
    /// If more than `u32::MAX` loci in one stratum share a single shape, or if the pooled base
    /// counters pass `u64`. Neither is reachable on a genome — the first needs 4.3 billion loci of
    /// one shape — so both are guards against a later edit rather than against an input. They are
    /// here because the release profile leaves `overflow-checks` off, so the alternative to a
    /// panic is a wrap: the stratum's commonest shape coming back holding no loci at all. The
    /// sibling path priced the same guard at 0.2 ns a site.
    pub fn add_locus(&mut self, shape: LocusShape, bases: BaseComparison) {
        let loci = self.loci_by_shape.entry(shape).or_insert(0);
        *loci = loci.checked_add(1).unwrap_or_else(|| {
            panic!("more than {} loci share one shape: {shape:?}", u32::MAX);
        });
        self.bases.add(bases.into());
    }

    /// Combine another shard's table into this one — element-wise integer addition, so it is
    /// associative and exact and a region-sharded walk merges to the table of the union.
    ///
    /// Destructured exhaustively, so a field added later is a compile error here rather than one
    /// that silently stops merging.
    ///
    /// # Panics
    ///
    /// If merging puts more than `u32::MAX` loci on one shape, or takes either pooled base
    /// counter past `u64`. Out of reach for any real genome, and kept for the reason
    /// [`Self::add_locus`]'s guards are.
    ///
    /// **Not atomic.** A table whose merge panicked has already absorbed the entries the loop
    /// reached. Nothing recovers from these panics, so nothing reads such a table; it is said here
    /// so that a later `catch_unwind` is not written on the assumption that it could.
    pub fn merge(&mut self, other: &Self) {
        let Self {
            loci_by_shape,
            bases,
        } = other;
        for (&shape, &loci) in loci_by_shape {
            let into = self.loci_by_shape.entry(shape).or_insert(0);
            *into = into.checked_add(loci).unwrap_or_else(|| {
                panic!(
                    "merging two shards puts more than {} loci on one shape: {shape:?}",
                    u32::MAX
                );
            });
        }
        self.bases.add(*bases);
    }

    /// Every entry — each distinct shape and how many loci had it — in shape order.
    ///
    /// **A `Vec` rather than an iterator**, because the search re-walks these once per candidate:
    /// a materialised list is a contiguous scan where re-walking the map is pointer-chasing
    /// through its nodes on every pass. It is also what keeps the storage decision private, which
    /// `&BTreeMap` would publish.
    #[must_use]
    pub fn entries(&self) -> Vec<StratumEntry> {
        self.loci_by_shape
            .iter()
            .map(|(&shape, &loci)| StratumEntry {
                shape,
                loci: u64::from(loci),
            })
            .collect()
    }

    /// How many distinct shapes this stratum holds, without building the list of them.
    #[must_use]
    pub fn entry_count(&self) -> usize {
        self.loci_by_shape.len()
    }

    /// How many loci this stratum holds — **loci, not entries**, and what `MIN_LOCI_TO_FIT` is
    /// read against. The two differ by however far the shapes collapsed: 2.3 loci an entry on
    /// HG002 at 300×, 25 on a whole tomato genome.
    ///
    /// Summed from the entries rather than kept beside them, so there is one place a locus is
    /// recorded and a merge cannot leave a count disagreeing with the table it describes. The sum
    /// cannot overflow: each entry is held below `u32::MAX` by the guard in [`Self::add_locus`],
    /// and a stratum cannot hold more than 646,645 entries — the number of shapes a locus capped
    /// at `MAX_LOCUS_READS` can take (`arch/parameter_prepass_ssr.md` §2.2) — which caps this at
    /// about 2.8 × 10¹⁵.
    #[must_use]
    pub fn loci(&self) -> u64 {
        self.loci_by_shape.values().copied().map(u64::from).sum()
    }

    /// The per-base substitution rate inside this stratum's tracts: mismatched bases over bases
    /// compared, and `None` where nothing was compared.
    ///
    /// **A division, and it is the maximum rather than a moment estimate** — a moment estimate
    /// being a number matched to an average rather than to the value that makes the data most
    /// likely. Each read's mismatches are binomial at one rate, so the log-likelihood over the
    /// pooled pair is maximised exactly at this ratio: a search over the rate finds the same
    /// number and has nothing to add (`spec/parameter_prepass_ssr.md` §4.1, and the test below
    /// that runs the search).
    ///
    /// **`None` rather than zero where nothing was compared**, because the two are different
    /// claims: a stratum of perfectly-matching reads has a measured rate of zero, and a stratum
    /// nothing was compared in has no rate at all. A zero would be borrowed from and fitted with
    /// as though it had been measured.
    ///
    /// Where a stratum holds reads of two different true rates, the pooled counters return their
    /// base-weighted mean, which is the right answer for a model carrying one rate.
    #[must_use]
    pub fn substitution_rate(&self) -> Option<ErrorRate> {
        if self.bases.compared == 0 {
            return None;
        }
        let rate = self.bases.mismatched as f64 / self.bases.compared as f64;
        // PANIC-FREE: the ratio is in `[0, 1]`, at every size and by two steps. `BaseComparison`
        // refuses a locus whose mismatched bases outnumber its compared ones, and both counters
        // are non-wrapping sums of such pairs — including through `merge`, whose argument is
        // another table holding the same invariant. And the conversion cannot break the
        // comparison: `u64 -> f64` is monotone, so `mismatched <= compared` survives it, and a
        // correctly-rounded quotient of `a <= b` cannot exceed one.
        Some(ErrorRate::try_new(rate).expect("mismatched bases never outnumber compared ones"))
    }

    /// Of the reads that differed from the reference tract length, the share that did so by
    /// something other than a whole number of motif copies — **the guard share**, which is what
    /// `GUARD_SHARE_LIMIT` bounds. Above it, what moved this stratum's reads is ordinary indel
    /// rather than repeat slippage, and a level fitted there is mostly mis-modelled indel however
    /// many loci stood behind it (`spec/parameter_prepass_ssr.md` §5).
    ///
    /// **The denominator is measured from the reference and the model means the allele**, which
    /// this table cannot know. A locus whose own allele differs from the reference contributes
    /// whole-repeat differences to this denominator and nothing to its numerator, so the reported
    /// share is diluted relative to the model's and never inflated: a stratum that crosses the
    /// limit on this number has crossed it on the model's too.
    ///
    /// Zero where no read differed from the reference at all — there is no share to report, and
    /// zero is the answer that leaves such a stratum unflagged, which is right: a stratum where
    /// nothing moved is not one this noise model fails on.
    ///
    /// Neither sum can overflow: both are bounded by the locus count times `MAX_LOCUS_READS`,
    /// which the entry-count and per-entry limits hold at about 3.3 × 10¹⁶.
    #[must_use]
    pub fn not_whole_repeat_share(&self) -> f64 {
        let mut off_reference = 0u64;
        let mut not_whole_repeat = 0u64;
        for (shape, &loci) in &self.loci_by_shape {
            let loci = u64::from(loci);
            off_reference += u64::from(shape.reads_off_reference()) * loci;
            not_whole_repeat += u64::from(shape.reads_not_whole_repeat()) * loci;
        }
        if off_reference == 0 {
            return 0.0;
        }
        not_whole_repeat as f64 / off_reference as f64
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

    // -----------------------------------------------------------------
    // The table of entries.
    // -----------------------------------------------------------------

    /// Bases compared and mismatched, for a fixture with no stake in the particular values.
    fn bases(bases_compared: u32, bases_mismatched: u32) -> BaseComparison {
        BaseComparison::try_new(bases_compared, bases_mismatched)
            .expect("mismatched within compared")
    }

    /// A table built from a sequence of loci, each a shape and its two base counts.
    fn table_of(loci: impl IntoIterator<Item = (LocusShape, BaseComparison)>) -> StratumTable {
        let mut table = StratumTable::default();
        for (shape, comparison) in loci {
            table.add_locus(shape, comparison);
        }
        table
    }

    /// One entry, for comparing against what a table returns.
    fn entry(shape: LocusShape, loci: u64) -> StratumEntry {
        StratumEntry { shape, loci }
    }

    /// **Two loci that looked alike are one entry with a count of two**, which is the collapse
    /// the whole design rests on: 1.73 million tomato loci made 70,305 entries, so all but 4 loci
    /// in 100 share a shape with another. A table keyed on anything unique to a locus would hold
    /// one entry each and the fit would be summing over loci rather than over shapes.
    #[test]
    fn two_loci_of_the_same_shape_make_one_entry_counted_twice() {
        let alike = shape_with_reads(&[(0, 4)], 0).expect("four reads");
        let different = shape_with_reads(&[(0, 3), (-1, 1)], 0).expect("four reads");

        let table = table_of([
            (alike, bases(400, 1)),
            (alike, bases(400, 2)),
            (different, bases(400, 0)),
        ]);

        assert_eq!(
            table.entries(),
            vec![entry(alike, 2), entry(different, 1)],
            "two entries, and the count says how many loci stood behind each"
        );
        assert_eq!(table.loci(), 3, "three loci, whatever they collapsed into");
        assert_eq!(table.entry_count(), 2);
    }

    /// **Loci and entries are different numbers**, and `MIN_LOCI_TO_FIT` is read against the
    /// first. A `loci()` that returned the entry count would let a stratum of 999 loci sharing
    /// one shape look like a stratum of one.
    #[test]
    fn the_locus_count_counts_loci_and_not_entries() {
        let one_shape = shape_with_reads(&[(0, 5)], 0).expect("five reads");
        let table = table_of([(one_shape, bases(500, 0)); 7]);

        assert_eq!(table.entry_count(), 1);
        assert_eq!(table.loci(), 7);
    }

    /// **The entries come back in an order the contents decide**, not the order the loci arrived
    /// in — which is what lets two runs over the same stratum return the same fitted level, since
    /// a fit is a floating-point sum over the entries and floating-point addition is not
    /// associative.
    ///
    /// **What this can fail on is the container, not the arithmetic.** A `BTreeMap`'s walk order
    /// is a property of the type, so no defect inside these methods could make the two orders
    /// differ; swapping the storage for a `HashMap` or a `Vec` is what it guards. The three loci
    /// carry *different* base counts so that the table equality beside it has power too — with
    /// one count repeated, a counter accumulating wrongly would agree with itself.
    #[test]
    fn the_entries_come_back_in_the_same_order_whatever_order_the_loci_arrived_in() {
        let loci = [
            (
                shape_with_reads(&[(0, 3)], 0).expect("three reads"),
                bases(300, 1),
            ),
            (
                shape_with_reads(&[(-4, 2), (0, 1)], 0).expect("three reads"),
                bases(280, 7),
            ),
            (
                shape_with_reads(&[(0, 2), (2, 1)], 1).expect("four reads"),
                bases(410, 0),
            ),
        ];

        let forwards = table_of(loci);
        let backwards = table_of(loci.into_iter().rev());

        assert_eq!(forwards.entries(), backwards.entries());
        assert_eq!(forwards, backwards);
    }

    /// **A stratum split across shards and merged is the stratum the whole walk would have
    /// built** — entry for entry and counter for counter, as an equality rather than a
    /// tolerance, which integer counts are what make possible. This is the property a
    /// region-sharded run rests on: cut the genome any way, merge in either order, get the same
    /// table.
    #[test]
    fn a_table_split_and_merged_equals_the_unsplit_one_in_either_order() {
        let loci = [
            (
                shape_with_reads(&[(0, 4)], 0).expect("four reads"),
                bases(400, 1),
            ),
            (
                shape_with_reads(&[(0, 3), (-1, 1)], 0).expect("four reads"),
                bases(400, 3),
            ),
            (
                shape_with_reads(&[(0, 4)], 0).expect("four reads"),
                bases(380, 0),
            ),
            (
                shape_with_reads(&[(1, 2)], 2).expect("four reads"),
                bases(410, 7),
            ),
            (
                shape_with_reads(&[(0, 3), (-1, 1)], 0).expect("four reads"),
                bases(400, 2),
            ),
        ];
        let whole = table_of(loci);

        for cut in 0..=loci.len() {
            let (left, right) = loci.split_at(cut);
            let mut forwards = table_of(left.iter().copied());
            forwards.merge(&table_of(right.iter().copied()));
            let mut backwards = table_of(right.iter().copied());
            backwards.merge(&table_of(left.iter().copied()));

            assert_eq!(
                forwards, whole,
                "split after {cut} loci, merged left to right"
            );
            assert_eq!(
                backwards, whole,
                "split after {cut} loci, merged right to left"
            );
        }
    }

    /// Three shards merged in each of the six orders give one table, and it is the table one
    /// walk would have built — which is what "associative and exact" means for a walk that does
    /// not promise which shard finishes first. **Compared against the whole table and not only
    /// across the orders**: any defect symmetric in merge order — dropping the base counters, say
    /// — is invisible to orders compared only with each other. Two of the shards share a shape,
    /// so the merge's add-into-an-existing-entry path is the one under test.
    #[test]
    fn three_shards_merge_to_the_same_table_in_every_order() {
        let loci = [
            (
                shape_with_reads(&[(0, 4)], 0).expect("four reads"),
                bases(400, 1),
            ),
            (
                shape_with_reads(&[(0, 3), (-1, 1)], 0).expect("four reads"),
                bases(390, 5),
            ),
            (
                shape_with_reads(&[(0, 4)], 0).expect("four reads"),
                bases(500, 0),
            ),
        ];
        let whole = table_of(loci);
        let shards: Vec<StratumTable> = loci.iter().map(|&locus| table_of([locus])).collect();

        for order in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            let mut merged = StratumTable::default();
            for shard in order {
                merged.merge(&shards[shard]);
            }
            assert_eq!(merged, whole, "merged in the order {order:?}");
        }
        assert_eq!(whole.loci(), 3);
        assert_eq!(whole.entry_count(), 2);
    }

    proptest::proptest! {
        /// **Any split of any list of loci merges back to the unsplit table.** The fixtures above
        /// fix the shard sizes and which shapes collide across the cut; this leaves both to the
        /// generator, over a shape space narrow enough that collisions are common, and includes
        /// the empty shard at either end.
        #[test]
        fn any_split_of_any_locus_list_merges_to_the_unsplit_table(
            loci in proptest::collection::vec(
                (0usize..OFFSET_BUCKETS, 1u32..4, 0u32..3, 0u32..500, 0u32..20),
                0..24usize,
            ),
            cut in 0usize..24,
        ) {
            let built: Vec<(LocusShape, BaseComparison)> = loci
                .iter()
                .map(|&(bucket, reads, guard, compared, mismatched)| {
                    let mut reads_by_bucket = [0u32; OFFSET_BUCKETS];
                    reads_by_bucket[bucket] = reads;
                    (
                        LocusShape::try_new(reads_by_bucket, guard).expect("inside the cap"),
                        BaseComparison::try_new(compared, mismatched.min(compared))
                            .expect("mismatched within compared"),
                    )
                })
                .collect();

            let whole = table_of(built.iter().copied());
            let cut = cut.min(built.len());
            let (left, right) = built.split_at(cut);

            let mut forwards = table_of(left.iter().copied());
            forwards.merge(&table_of(right.iter().copied()));
            let mut backwards = table_of(right.iter().copied());
            backwards.merge(&table_of(left.iter().copied()));

            proptest::prop_assert_eq!(&forwards, &whole);
            proptest::prop_assert_eq!(&backwards, &whole);
            proptest::prop_assert_eq!(forwards.loci(), built.len() as u64);
        }
    }

    /// **The overflow guard is what stands between a wrapped count and a stratum whose commonest
    /// shape reports no loci at all**, and the release profile leaves `overflow-checks` off, so
    /// nothing else would. Reached by merging a table with a copy of itself, which doubles every
    /// count: thirty-two doublings pass `u32::MAX`.
    #[test]
    #[should_panic(expected = "merging two shards puts more than")]
    fn merging_past_a_u32_of_loci_on_one_shape_panics_rather_than_wrapping() {
        let shape = shape_with_reads(&[(0, 1)], 0).expect("one read");
        let mut table = table_of([(shape, bases(0, 0))]);

        for _ in 0..32 {
            let copy = table.clone();
            table.merge(&copy);
        }
    }

    /// The same limit reached one locus at a time, which is the guard on the other entry point.
    /// `2^0 + 2^1 + … + 2^31` is exactly `u32::MAX` loci on one shape, so the next `add_locus` is
    /// the one that would wrap.
    #[test]
    #[should_panic(expected = "loci share one shape")]
    fn adding_a_locus_past_a_u32_of_them_on_one_shape_panics_rather_than_wrapping() {
        let shape = shape_with_reads(&[(0, 1)], 0).expect("one read");
        let doubled = |power: u32| {
            let mut table = table_of([(shape, bases(0, 0))]);
            for _ in 0..power {
                let copy = table.clone();
                table.merge(&copy);
            }
            table
        };

        let mut table = StratumTable::default();
        for power in 0..32 {
            table.merge(&doubled(power));
        }
        assert_eq!(
            table.loci(),
            u64::from(u32::MAX),
            "one locus short of the guard"
        );

        table.add_locus(shape, bases(0, 0));
    }

    /// **The base counters have their own guard**, and a wrapped one is worse than a wrapped
    /// locus count: the stratum then reports a substitution rate computed from a remainder. Four
    /// distinct shapes at `u32::MAX` compared bases each, doubled 31 times and merged, so the
    /// base counters reach `u64` before any shape's locus count reaches `u32`.
    #[test]
    #[should_panic(expected = "bases compared in one stratum passed")]
    fn merging_past_a_u64_of_compared_bases_panics_rather_than_wrapping() {
        let half = |first: LocusShape, second: LocusShape| {
            let mut table = table_of([
                (first, bases(u32::MAX, u32::MAX)),
                (second, bases(u32::MAX, u32::MAX)),
            ]);
            for _ in 0..31 {
                let copy = table.clone();
                table.merge(&copy);
            }
            table
        };

        let mut left = half(
            shape_with_reads(&[(0, 1)], 0).expect("one read"),
            shape_with_reads(&[(1, 1)], 0).expect("one read"),
        );
        let right = half(
            shape_with_reads(&[(2, 1)], 0).expect("one read"),
            shape_with_reads(&[(3, 1)], 0).expect("one read"),
        );

        left.merge(&right);
    }

    /// **A stratum where nothing was compared has no substitution rate**, and that is a
    /// different claim from a rate of zero: zero says every base read matched, `None` says no
    /// base was read. A zero would be borrowed from and fitted with as though it had been
    /// measured.
    #[test]
    fn a_table_with_no_bases_compared_has_no_substitution_rate() {
        assert_eq!(StratumTable::default().substitution_rate(), None);

        let no_bases = table_of([(
            shape_with_reads(&[(0, 4)], 0).expect("four reads"),
            bases(0, 0),
        )]);
        assert_eq!(no_bases.substitution_rate(), None);

        let matched = table_of([(
            shape_with_reads(&[(0, 4)], 0).expect("four reads"),
            bases(400, 0),
        )]);
        assert_eq!(
            matched.substitution_rate().map(ErrorRate::get),
            Some(0.0),
            "every base matching is a measured rate of zero, not an absent one"
        );
    }

    /// **A stratum every one of whose compared bases mismatched reports a rate of one**, which is
    /// the upper edge the `expect` in `substitution_rate` claims to be safe at. It is safe
    /// because `ErrorRate` accepts both endpoints — a fact in another module that nothing here
    /// would otherwise pin, and whose narrowing would land as a panic on the noisiest stratum in
    /// a run.
    #[test]
    fn a_stratum_whose_every_compared_base_mismatched_reports_a_rate_of_one() {
        let table = table_of([(
            shape_with_reads(&[(0, 4)], 0).expect("four reads"),
            bases(400, 400),
        )]);

        assert_eq!(
            table.substitution_rate().map(ErrorRate::get),
            Some(1.0),
            "every base mismatching is a measured rate of one, not a panic"
        );
    }

    /// **The closed form is the maximum, not merely a ratio**, which is the claim that lets this
    /// path fit its fourth number by division while the other three are searched.
    ///
    /// Each read's mismatches are binomial at one rate, so a stratum's log-likelihood is
    /// `m·ln ε + (c−m)·ln(1−ε)`. This walks that curve at 99,999 interior rates, a hundred-
    /// thousandth apart, and asks two things of the rate the **table** returned: that the grid's
    /// best sits within one step of it, and that no grid point scores above it. Four
    /// (compared, mismatched) pairs, spanning a rate of 1 in 400 to one of 63 in 64, so the
    /// answer is not one arithmetic coincidence.
    ///
    /// **The grid is compared against what the code returns, not against a literal.** An earlier
    /// version pinned the table to `0.0030` within `1e-12` first, which is ten million times
    /// tighter than a grid step — so the grid assertions could never be the ones to fail, and the
    /// oracle proved nothing about the code.
    #[test]
    fn the_substitution_rate_is_the_maximum_at_every_mismatch_count_it_was_tried_at() {
        let shape = shape_with_reads(&[(0, 4)], 0).expect("four reads");

        for (compared, mismatched) in [(10_000u32, 30u32), (1_000, 500), (400, 1), (64, 63)] {
            let table = table_of([(shape, bases(compared, mismatched))]);
            let closed_form = table
                .substitution_rate()
                .expect("bases were compared")
                .get();

            let matched = f64::from(compared - mismatched);
            let missed = f64::from(mismatched);
            let score = |rate: f64| missed * rate.ln() + matched * (1.0 - rate).ln();

            let steps = 100_000;
            let (best_rate, best_score) = (1..steps)
                .map(|step| {
                    let rate = f64::from(step) / f64::from(steps);
                    (rate, score(rate))
                })
                .fold((f64::NAN, f64::NEG_INFINITY), |best, here| {
                    if here.1 > best.1 { here } else { best }
                });

            assert!(
                (best_rate - closed_form).abs() <= 1.0 / f64::from(steps),
                "at {mismatched} mismatches in {compared} bases the searched maximum {best_rate} \
                 is more than one grid step from the closed form {closed_form}"
            );
            assert!(
                score(closed_form) >= best_score,
                "at {mismatched} mismatches in {compared} bases a grid point scored {best_score}, \
                 above the closed form's {}",
                score(closed_form)
            );
        }
    }

    /// The two rates the grid above cannot reach, because there the maximum sits **on** the edge
    /// of `[0, 1]` rather than inside it: every base matching puts it at zero, every base
    /// mismatching at one. Stated as the answer rather than as a search, since the score is
    /// undefined at both endpoints — `0 · ln 0` is not a number.
    #[test]
    fn the_substitution_rate_sits_on_the_boundary_when_every_base_agrees_or_every_base_differs() {
        let shape = shape_with_reads(&[(0, 4)], 0).expect("four reads");

        let all_matched = table_of([(shape, bases(400, 0))]);
        assert_eq!(
            all_matched.substitution_rate().map(ErrorRate::get),
            Some(0.0)
        );

        let all_mismatched = table_of([(shape, bases(400, 400))]);
        assert_eq!(
            all_mismatched.substitution_rate().map(ErrorRate::get),
            Some(1.0)
        );
    }

    /// Two shards' base counts pool into one rate, and the rate a stratum reports is the
    /// base-weighted mean of what its reads actually carried — which is the right answer for a
    /// model carrying one rate per stratum. The two loci are far apart, 1 mismatch in 1,000
    /// against 40, so a locus-weighted mean (0.0245) and a base-weighted one (0.0049) cannot be
    /// confused.
    #[test]
    fn merged_base_counts_give_the_base_weighted_rate() {
        let shape = shape_with_reads(&[(0, 4)], 0).expect("four reads");
        let mut clean = table_of([(shape, bases(9_000, 9))]);
        let noisy = table_of([(shape, bases(1_000, 40))]);

        clean.merge(&noisy);

        let rate = clean.substitution_rate().expect("bases compared").get();
        assert!(
            (rate - 49.0 / 10_000.0).abs() < 1e-12,
            "49 mismatches in 10,000 bases is 0.0049, not {rate}"
        );
    }

    /// **The guard share is over the reads that moved, not over every read**, which is what
    /// makes it a statement about how this stratum's reads moved rather than about how much they
    /// moved at all. Two loci, twenty reads, four off the reference length of which one is
    /// non-whole-repeat: the share is one in four, where a denominator of every read would give
    /// one in twenty.
    #[test]
    fn the_guard_share_is_measured_over_the_reads_that_left_the_reference_length() {
        let table = table_of([
            (
                shape_with_reads(&[(0, 7), (-1, 3)], 0).expect("ten reads"),
                bases(1_000, 3),
            ),
            (
                shape_with_reads(&[(0, 9)], 1).expect("ten reads"),
                bases(1_000, 3),
            ),
        ]);

        assert!(
            (table.not_whole_repeat_share() - 0.25).abs() < 1e-12,
            "one guard read among four off-reference reads is 0.25, not {}",
            table.not_whole_repeat_share()
        );
    }

    /// **A locus counted by two is a locus counted twice in the share**, so a stratum whose
    /// commonest shape is interrupted reports that, rather than one entry's worth of it.
    #[test]
    fn the_guard_share_weighs_each_entry_by_how_many_loci_had_it() {
        let interrupted = shape_with_reads(&[(0, 8)], 2).expect("ten reads");
        let clean = shape_with_reads(&[(0, 8), (-1, 2)], 0).expect("ten reads");

        let one_each = table_of([(interrupted, bases(1_000, 1)), (clean, bases(1_000, 1))]);
        assert!(
            (one_each.not_whole_repeat_share() - 0.5).abs() < 1e-12,
            "two guard reads among four off-reference: {}",
            one_each.not_whole_repeat_share()
        );

        let mut interrupted_thrice = one_each.clone();
        for _ in 0..2 {
            interrupted_thrice.add_locus(interrupted, bases(1_000, 1));
        }
        assert!(
            (interrupted_thrice.not_whole_repeat_share() - 0.75).abs() < 1e-12,
            "six guard reads among eight off-reference: {}",
            interrupted_thrice.not_whole_repeat_share()
        );
    }

    /// **A real non-reference allele lands in the denominator and not the numerator**, so the
    /// share this table reports is diluted relative to the model's and never inflated — which is
    /// what makes a stratum crossing `GUARD_SHARE_LIMIT` here certain to have crossed it on the
    /// model's number too. Both tables hold the same two guard reads; the second adds a locus
    /// whose ten reads all sit a whole repeat below the reference, and its share falls from 2 in
    /// 2 to 2 in 12.
    #[test]
    fn a_non_reference_allele_dilutes_the_guard_share_rather_than_inflating_it() {
        let interrupted = shape_with_reads(&[(0, 8)], 2).expect("ten reads");
        let shifted_allele = shape_with_reads(&[(-1, 10)], 0).expect("ten reads, a repeat short");

        let guard_only = table_of([(interrupted, bases(1_000, 1))]);
        let with_allele = table_of([
            (interrupted, bases(1_000, 1)),
            (shifted_allele, bases(1_000, 1)),
        ]);

        assert!(
            (guard_only.not_whole_repeat_share() - 1.0).abs() < 1e-12,
            "every off-reference read here is a guard read"
        );
        assert!(
            (with_allele.not_whole_repeat_share() - 2.0 / 12.0).abs() < 1e-12,
            "the allele's ten reads join the denominator alone: {}",
            with_allele.not_whole_repeat_share()
        );
    }

    /// A stratum whose reads *did* leave the reference length, every one of them by a whole
    /// number of copies, reports a share of zero — which is the answer that separates "nothing
    /// moved" from "everything that moved moved by whole repeats". A numerator borrowing from
    /// the denominator would report something else.
    #[test]
    fn a_stratum_whose_off_reference_reads_are_all_whole_repeat_reports_a_zero_share() {
        let table = table_of([(
            shape_with_reads(&[(0, 4), (-1, 3), (2, 3)], 0).expect("ten reads"),
            bases(1_000, 2),
        )]);

        assert_eq!(table.not_whole_repeat_share(), 0.0);
    }

    /// A stratum where every read sat at the reference length has no share to report, and zero
    /// is the answer that leaves it unflagged — which is right, since nothing moved for this
    /// noise model to describe wrongly.
    #[test]
    fn a_stratum_where_no_read_left_the_reference_length_reports_no_guard_share() {
        let table = table_of([(
            shape_with_reads(&[(0, 10)], 0).expect("ten reads"),
            bases(1_000, 2),
        )]);

        assert_eq!(table.not_whole_repeat_share(), 0.0);
        assert_eq!(StratumTable::default().not_whole_repeat_share(), 0.0);
    }

    /// **More mismatched bases than compared ones is a swapped pair, not a noisy locus**, and
    /// refusing it is what makes the fitted rate a probability by construction.
    #[test]
    fn a_locus_cannot_offer_more_mismatched_bases_than_it_compared() {
        assert_eq!(
            BaseComparison::try_new(3, 400),
            Err(DomainError::SsrBaseComparison {
                bases_compared: 3,
                bases_mismatched: 400
            }),
            "the two counts swapped"
        );
        assert!(
            BaseComparison::try_new(400, 400).is_ok(),
            "every base mismatching is a real, if implausible, locus"
        );
    }

    /// **The refusal bites at the first pair it should**, one mismatch past the bases compared,
    /// and not one pair later. A guard off by one there admits a locus whose rate exceeds one,
    /// and the run then dies inside `substitution_rate` — at the `expect` whose comment says it
    /// cannot.
    #[test]
    fn a_base_comparison_is_refused_at_one_mismatch_past_the_bases_compared() {
        assert_eq!(
            BaseComparison::try_new(400, 401),
            Err(DomainError::SsrBaseComparison {
                bases_compared: 400,
                bases_mismatched: 401
            }),
            "one more mismatch than compared is the first pair to refuse"
        );
        assert_eq!(
            BaseComparison::try_new(0, 1),
            Err(DomainError::SsrBaseComparison {
                bases_compared: 0,
                bases_mismatched: 1
            }),
            "the same boundary where nothing was compared"
        );
        assert!(
            BaseComparison::try_new(0, 0).is_ok(),
            "a locus that compared nothing mismatched nothing"
        );
    }

    /// Each count comes back as it was given — so a transposed accessor pair fails *here*, at
    /// the two-line accessor, rather than three call sites downstream as a substitution rate of
    /// 333.
    #[test]
    fn a_base_comparison_reports_each_of_its_two_counts_as_it_was_given_them() {
        let pair = BaseComparison::try_new(400, 7).expect("seven within four hundred");

        assert_eq!(pair.bases_compared(), 400);
        assert_eq!(pair.bases_mismatched(), 7);
    }

    /// An empty table is a stratum nothing was entered into, and every accessor says so rather
    /// than dividing by zero or reporting a number nothing stands behind.
    #[test]
    fn an_empty_table_reports_no_loci_no_entries_and_no_rate() {
        let empty = StratumTable::default();

        assert_eq!(empty.loci(), 0);
        assert_eq!(empty.entry_count(), 0);
        assert!(empty.entries().is_empty());
        assert_eq!(empty.substitution_rate(), None);
        assert_eq!(empty.not_whole_repeat_share(), 0.0);
    }
}
