//! The cell table: a tally of what the sites looked like.
//!
//! Each covered position reduces to a small key — how many reads covered it, how many
//! of those showed something other than the reference, and, when a sample has more
//! than one read group and few alternative reads, which library those reads came
//! from. The table counts how many positions showed each key. Eight hundred million
//! positions become a few hundred counters, and the fits lose nothing, because a site
//! enters the likelihood only through that key
//! (`spec/parameter_prepass_generic.md` §4): the positions that all looked alike are
//! scored once and multiplied by how many there were.
//!
//! **The library attribution is the part that is easy to drop and cannot be.** With
//! it forgotten, a key of total depth and total alternative count sees only the
//! share-weighted mean error rate and nothing else about the individual libraries —
//! the likelihood is exactly flat along every combination holding that mean fixed, so
//! no amount of genome separates them. Keeping which library each of the first few
//! alternative reads came from is what breaks that flatness
//! (`arch/parameter_prepass_generic.md` §2.2).
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §2.2. The table itself
//! and the fold over it land across Milestone B; the depth ladder it is binned by is
//! [`super::depth_bins`], from Milestone A.

use crate::ng::types::ReadGroupId;

use smallvec::SmallVec;

use super::depth_bins::DepthBin;

/// How many alternative reads keep their library attribution. Above this the site
/// pools them.
///
/// **A precision choice, not a correctness one**, which is what changed when the
/// scoring rule became a likelihood: the fit is exactly unbiased at any value of this
/// bound, so what it trades is cells against how sharply two libraries can be told
/// apart. Four is the measured default — at three reads a bound of two is equally
/// unbiased on 28% fewer cells, and neither loses measurable precision against scoring
/// every read against its own library
/// (`spec/parameter_prepass_generic.md` §1, §12.8).
pub const MAX_ATTRIBUTED_ALT_READS: u32 = 4;

// A site's attributed counts are held one byte each, since none of them can exceed the
// bound above. Checked here rather than by test, because the cast that relies on it is
// infallible only while this holds.
const _: () = assert!(
    MAX_ATTRIBUTED_ALT_READS <= u8::MAX as u32,
    "an attributed alternative count is stored in one byte"
);

/// One site's evidence, reduced: how many reads covered it, and how many of those
/// showed something other than the reference.
///
/// Named for the two numbers it holds rather than for where they came from, so that a
/// call site reading `-> DepthAndAltReads` needs nothing else to know what it got.
///
/// **`alt_reads <= depth` always, and the constructor is what makes that true.** It is
/// the invariant the whole depth-per-cell argument rests on: every site in a cell has
/// at least as many reads as it has alternative ones, so the cell's mean depth does
/// too, and no cell is ever charged a negative count of reference reads
/// (`arch/parameter_prepass_generic.md` §2.2). A transposed pair would satisfy
/// nothing downstream and would be scored as a plausible site.
///
/// **It asserts rather than returning a `Result`, unlike the constrained scalars in
/// `types.rs`.** Those guard values that arrive from outside — a BED, a flag, a header
/// — where a caller can report and carry on. These two are counted from a locus's own
/// reads by [`super::depth_and_alt_reads`], so a violation is this crate's arithmetic
/// being broken rather than a data condition, there is no second thing a caller could
/// do about it, and the check sits on the path walked once per covered position.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct DepthAndAltReads {
    depth: u32,
    alt_reads: u32,
}

impl DepthAndAltReads {
    /// The only constructor.
    ///
    /// # Panics
    ///
    /// If `alt_reads` exceeds `depth` — see the type's doc for why that is a panic.
    #[must_use]
    pub fn new(depth: u32, alt_reads: u32) -> Self {
        assert!(
            alt_reads <= depth,
            "a site of depth {depth} cannot show {alt_reads} alternative reads — the \
             pair is transposed or was counted from different read sets"
        );
        Self { depth, alt_reads }
    }

    /// How many reads covered the site and witnessed something there.
    #[inline]
    #[must_use]
    pub fn depth(self) -> u32 {
        self.depth
    }

    /// How many of those reads showed something other than the reference.
    #[inline]
    #[must_use]
    pub fn alt_reads(self) -> u32 {
        self.alt_reads
    }
}

/// The tally's key: what a site looked like, at the grain the fits read it.
///
/// **One entry per site** — this is a tally, and sites with identical keys are counted
/// together.
///
/// **Two arms, and which one a site takes is decided here.** The alternative reads
/// keep the library each came from while there are at most
/// [`MAX_ATTRIBUTED_ALT_READS`] of them; above that everything pools, and those are
/// the cells where the genotype is certain and the error rate is not being estimated
/// from them anyway. At one read group both arms say the same thing, so every key is
/// pooled and 1,550 of the 1,707 samples in the archive survey are keyed exactly as
/// they would be with no multi-library machinery at all.
///
/// **The attributed arm keeps a depth bin and not an exact depth, and that is measured
/// rather than assumed.** Under the adopted ladder the fit's asymptotic bias stays at
/// 0.054 rungs of the error-rate ladder and 0.3% in each genotype frequency, against
/// exactly zero unbinned (research note §4.3).
///
/// **A struct with private fields rather than the public two-variant enum the
/// architecture sketches** (`arch/parameter_prepass_generic.md` §2.2), because two of
/// this key's properties are invariants whose violation is silent and neither can be
/// expressed on an enum — a variant's fields are exactly as public as the enum itself:
///
/// - **The listing is canonical.** Two sites whose read groups were listed in
///   different orders, or one of which named a group that contributed no alternative
///   read, are the same site and must reach the same cell. Unsorted keys would split
///   one cell into several — each scored correctly, so nothing downstream would show
///   it, while `cells()`'s count and the region-sharding identity of
///   `spec/parameter_prepass_generic.md` §12.6 would both move.
/// - **Which arm a site takes follows from its counts**, so it is not a choice a call
///   site gets to make differently from its neighbour.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct SiteKey {
    /// Which depths this site is counted with.
    depth_bin: DepthBin,
    /// How many reads showed something other than the reference, whether or not the
    /// attribution below survived.
    alt_reads: u32,
    /// Which library each alternative read came from, ascending by read group and
    /// carrying no zero entries. **Empty exactly when the arm is pooled**, which the
    /// constructors below are what guarantee.
    alt_by_group: SmallVec<[(ReadGroupId, u8); 2]>,
}

impl SiteKey {
    /// The key for a site whose alternative reads are not attributed to a library:
    /// every key of the read-group histogram, where the whole entry belongs to one
    /// group already, and every key of a single-library sample.
    #[must_use]
    pub fn pooled(depth_bin: DepthBin, alt_reads: u32) -> Self {
        Self {
            depth_bin,
            alt_reads,
            alt_by_group: SmallVec::new(),
        }
    }

    /// The key for a site of a multi-library sample, keeping which library each
    /// alternative read came from where that is worth keeping.
    ///
    /// `alt_by_group` is one entry per read group that showed an alternative read at
    /// this site, in any order; the key it produces is canonical. Three things pool
    /// instead, and each is the same information as an attributed key rather than less
    /// of it:
    ///
    /// - **No alternative read at all.** There is no attribution to keep, and this is
    ///   the cell the overwhelming majority of the genome lands in — it belongs in the
    ///   dense table with its neighbours, not in the sparse map.
    /// - **More than [`MAX_ATTRIBUTED_ALT_READS`] of them.**
    /// - **A group that showed none.** Its absence from the listing says the same
    ///   thing its zero entry would, so carrying the zero would give one site two
    ///   possible keys.
    ///
    /// # Panics
    ///
    /// If a read group appears twice — that is one site counted from two overlapping
    /// read sets, and summing the duplicates would hide it.
    #[must_use]
    pub fn attributing(depth_bin: DepthBin, alt_by_group: &[(ReadGroupId, u32)]) -> Self {
        let alt_reads = alt_by_group.iter().fold(0u32, |total, &(group, alt)| {
            total.checked_add(alt).unwrap_or_else(|| {
                panic!(
                    "read group {} pushed a site's alternative count past u32",
                    group.get()
                )
            })
        });

        if alt_reads == 0 || alt_reads > MAX_ATTRIBUTED_ALT_READS {
            return Self::pooled(depth_bin, alt_reads);
        }

        let mut attributed: SmallVec<[(ReadGroupId, u8); 2]> = alt_by_group
            .iter()
            .filter(|&&(_, alt)| alt > 0)
            // PANIC-FREE: the counts sum to at most `MAX_ATTRIBUTED_ALT_READS` here,
            // which the const assertion above pins inside one byte.
            .map(|&(group, alt)| (group, alt as u8))
            .collect();
        attributed.sort_unstable_by_key(|&(group, _)| group);
        assert!(
            attributed.windows(2).all(|pair| pair[0].0 < pair[1].0),
            "a site listed the same read group twice: {attributed:?}"
        );

        Self {
            depth_bin,
            alt_reads,
            alt_by_group: attributed,
        }
    }

    /// Which depths this site is counted with.
    #[inline]
    #[must_use]
    pub fn depth_bin(&self) -> DepthBin {
        self.depth_bin
    }

    /// How many reads showed something other than the reference — the total, on
    /// either arm.
    #[inline]
    #[must_use]
    pub fn alt_reads(&self) -> u32 {
        self.alt_reads
    }

    /// Which library each alternative read came from, ascending by read group, or
    /// `None` where the site pooled them.
    ///
    /// `None` is not "this sample has one library": it is "this key does not say", and
    /// the score for such a cell sums over the splits it forgot rather than inventing
    /// one (`arch/parameter_prepass_generic.md` §5.1).
    #[inline]
    #[must_use]
    pub fn attribution(&self) -> Option<&[(ReadGroupId, u8)]> {
        (!self.alt_by_group.is_empty()).then_some(&self.alt_by_group)
    }
}

/// What a cell may be counted in.
///
/// **Two implementors, `u32` and `u64`, and the trait is sealed so there is no
/// third.** It exists to make the widening at the whole-sample fold explicit in the
/// type rather than to invite arithmetic over "any integer": a window accumulates in
/// `u32` and the fold over a sample's windows produces `u64`, because the two hold
/// quantities four orders of magnitude apart
/// (`arch/parameter_prepass_generic.md` §2.2).
///
/// **Addition is checked, and that is the point of having the trait at all.** This
/// repo's release profile leaves `overflow-checks` off, so a plain `+=` on a counter
/// that ran past its width would wrap to a small number and be scored as a rare cell —
/// silently, in the one module whose whole difficulty is that its wrong numbers have
/// no symptom. A panic naming the cell is the alternative, and it costs one predictable
/// branch per site.
pub trait CellCounter: Copy + Into<u64> + From<u32> + std::fmt::Debug + sealed::Sealed {
    /// A cell nothing has landed in.
    const ZERO: Self;

    /// Add, or `None` if the sum would not fit. The caller supplies the message,
    /// because what a reader needs is which cell overflowed and not which integer type.
    #[must_use]
    fn try_add(self, other: Self) -> Option<Self>;
}

impl CellCounter for u32 {
    const ZERO: Self = 0;

    #[inline]
    fn try_add(self, other: Self) -> Option<Self> {
        u32::checked_add(self, other)
    }
}

impl CellCounter for u64 {
    const ZERO: Self = 0;

    #[inline]
    fn try_add(self, other: Self) -> Option<Self> {
        u64::checked_add(self, other)
    }
}

mod sealed {
    /// Closes [`super::CellCounter`] to the two widths the design names. A third
    /// implementor would be a width nothing in `spec/parameter_prepass_generic.md` §9
    /// prices and that the fold has no rule for.
    pub trait Sealed {}
    impl Sealed for u32 {}
    impl Sealed for u64 {}
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(id: u32) -> ReadGroupId {
        ReadGroupId(id)
    }

    /// A depth and an alternative count are two `u32`s in the same order of magnitude,
    /// so the pair is the one place that can tell them apart — and it does, by
    /// refusing the transposition rather than storing it.
    #[test]
    fn a_site_cannot_show_more_alternative_reads_than_it_has_reads() {
        let site = DepthAndAltReads::new(12, 1);
        assert_eq!(site.depth(), 12);
        assert_eq!(site.alt_reads(), 1);

        // The boundary is legal: every read showed the alternative allele.
        let homozygous = DepthAndAltReads::new(12, 12);
        assert_eq!(homozygous.depth(), homozygous.alt_reads());
    }

    #[test]
    #[should_panic(expected = "cannot show 12 alternative reads")]
    fn a_transposed_depth_and_alternative_count_is_refused() {
        let _ = DepthAndAltReads::new(1, 12);
    }

    /// **The canonical-listing property, which is what keeps one site out of two
    /// cells.** `count_by_read_group` walks a locus's reads and emits a group's entry
    /// when it first meets one, so the order it produces follows the reads rather than
    /// the read-group table — two sites with identical evidence can arrive listed
    /// differently, and unsorted they would be counted apart.
    #[test]
    fn two_sites_listing_their_read_groups_in_either_order_reach_the_same_cell() {
        let bin = DepthBin(12);
        let forwards = SiteKey::attributing(bin, &[(group(0), 1), (group(7), 2)]);
        let backwards = SiteKey::attributing(bin, &[(group(7), 2), (group(0), 1)]);

        assert_eq!(forwards, backwards);
        assert_eq!(
            forwards.attribution(),
            Some(&[(group(0), 1u8), (group(7), 2u8)][..]),
            "ascending by read group"
        );
        assert_eq!(forwards.alt_reads(), 3);
    }

    /// A read group that covered the site and showed nothing says the same thing by
    /// being absent, so carrying its zero would give one site two possible keys.
    #[test]
    fn a_read_group_showing_no_alternative_read_is_dropped_rather_than_listed_as_zero() {
        let bin = DepthBin(12);
        let listed = SiteKey::attributing(bin, &[(group(0), 2), (group(7), 0)]);
        let omitted = SiteKey::attributing(bin, &[(group(0), 2)]);

        assert_eq!(listed, omitted);
        assert_eq!(listed.attribution(), Some(&[(group(0), 2u8)][..]));
    }

    /// **The two ends of the attributed arm both pool, and for different reasons.**
    /// With no alternative read there is nothing to attribute, and that cell holds the
    /// overwhelming majority of the genome — it belongs in the dense table. Above the
    /// bound the genotype is certain and the error rate is not estimated from those
    /// cells anyway.
    #[test]
    fn a_site_with_nothing_to_attribute_or_too_much_pools() {
        let bin = DepthBin(19);

        let nothing = SiteKey::attributing(bin, &[(group(0), 0), (group(7), 0)]);
        assert_eq!(nothing, SiteKey::pooled(bin, 0));
        assert_eq!(nothing.attribution(), None);

        let at_the_bound = SiteKey::attributing(bin, &[(group(0), 2), (group(7), 2)]);
        assert_eq!(at_the_bound.alt_reads(), MAX_ATTRIBUTED_ALT_READS);
        assert!(
            at_the_bound.attribution().is_some(),
            "the bound is inclusive, so four alternative reads keep their libraries"
        );

        let over_the_bound = SiteKey::attributing(bin, &[(group(0), 2), (group(7), 3)]);
        assert_eq!(over_the_bound, SiteKey::pooled(bin, 5));
        assert_eq!(over_the_bound.attribution(), None);
        assert_eq!(
            over_the_bound.alt_reads(),
            5,
            "pooling keeps the count it could not attribute"
        );
    }

    /// At one read group the attributed arm carries nothing the pooled key does not,
    /// and the design says so in as many words — but a lone library's site still takes
    /// the attributed arm when it shows a few alternative reads, because the *sample*
    /// is what decides whether attribution is worth keeping, and that decision belongs
    /// to the accumulator (Milestone C), not to this key.
    #[test]
    fn a_lone_read_group_still_produces_a_distinct_key_from_the_pooled_one() {
        let bin = DepthBin(12);
        let attributed = SiteKey::attributing(bin, &[(group(3), 1)]);

        assert_ne!(attributed, SiteKey::pooled(bin, 1));
        assert_eq!(attributed.alt_reads(), SiteKey::pooled(bin, 1).alt_reads());
    }

    #[test]
    #[should_panic(expected = "listed the same read group twice")]
    fn a_site_listing_one_read_group_twice_is_refused() {
        let _ = SiteKey::attributing(DepthBin(12), &[(group(4), 1), (group(4), 1)]);
    }

    /// The key sorts bin-major, then by alternative count, then by attribution — which
    /// is what makes a `BTreeMap` of these iterate in an order that does not vary
    /// between runs, and every fit is a floating-point sum over cells, which is not
    /// associative.
    #[test]
    fn keys_order_by_depth_bin_then_alternative_count() {
        let mut keys = [
            SiteKey::pooled(DepthBin(12), 3),
            SiteKey::attributing(DepthBin(9), &[(group(1), 1)]),
            SiteKey::pooled(DepthBin(9), 1),
            SiteKey::pooled(DepthBin(9), 0),
        ];
        keys.sort();

        let shape: Vec<(u16, u32, bool)> = keys
            .iter()
            .map(|key| {
                (
                    key.depth_bin().get(),
                    key.alt_reads(),
                    key.attribution().is_some(),
                )
            })
            .collect();
        assert_eq!(
            shape,
            vec![(9, 0, false), (9, 1, false), (9, 1, true), (12, 3, false)]
        );
    }

    /// The bound is what the research note measured, so it is stated as a literal
    /// rather than left to be read off a cell count.
    #[test]
    fn the_attribution_bound_is_the_measured_one() {
        assert_eq!(MAX_ATTRIBUTED_ALT_READS, 4);
    }

    /// **A counter that ran past its width must panic, not wrap.** The release profile
    /// this repo builds leaves `overflow-checks` off, so a wrapped depth sum would come
    /// back as a small number and be scored as a shallow cell — which is the failure
    /// the per-cell depth sum exists to prevent, arriving by another route.
    #[test]
    fn a_counter_reports_the_sum_it_cannot_hold_rather_than_wrapping() {
        assert_eq!(CellCounter::try_add(1u32, 2u32), Some(3u32));
        assert_eq!(CellCounter::try_add(u32::MAX, 1u32), None);

        // The width the whole-sample fold widens to holds what `u32` cannot: a human
        // genome's depth sum reaches 3.1 × 10¹¹ against a `u32` ceiling of 4.29 × 10⁹.
        let human_depth_sum = 310_000_000_000u64;
        assert!(human_depth_sum > u64::from(u32::MAX));
        assert_eq!(
            CellCounter::try_add(human_depth_sum, 1),
            Some(310_000_000_001)
        );
        assert_eq!(CellCounter::try_add(u64::MAX, 1u64), None);
    }

    /// Both widths start empty and widen upward without loss, which is what the fold
    /// of Milestone B4 rests on.
    #[test]
    fn both_counter_widths_start_at_zero_and_widen_to_u64() {
        assert_eq!(<u32 as CellCounter>::ZERO, 0);
        assert_eq!(<u64 as CellCounter>::ZERO, 0);
        assert_eq!(Into::<u64>::into(u32::MAX), 4_294_967_295u64);
        assert_eq!(<u64 as From<u32>>::from(7), 7u64);
    }
}
