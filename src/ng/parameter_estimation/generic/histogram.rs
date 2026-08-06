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
//! **Two words that are not the same word.** A **locus** is one entry in the tally —
//! one key, counted once. A **position** is one reference base. They differ because a
//! generic locus can be widened to an indel's reference span, which makes it one locus
//! and several positions, and the two counters below report them separately for exactly
//! that reason (`spec/parameter_prepass_generic.md` §4). "Site" is used here for a locus,
//! following the architecture's own names — [`SiteKey`],
//! [`DepthAltHistogram::add_site`] — and never for a position.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §2.2. The table itself
//! and the fold over it land across Milestone B; the depth ladder it is binned by is
//! [`super::depth_bins`], from Milestone A.

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::ng::types::{Bp, Ploidy, ReadGroupId};

use smallvec::SmallVec;

use super::depth_bins::{DepthBin, DepthBinEdges};

/// How many alternative reads keep their library attribution. Above this the site
/// pools them.
///
/// **A precision choice, not a correctness one**, which is what changed when the
/// scoring rule became a likelihood: the fit is unbiased at any bound that keeps some
/// attribution, so what the bound trades is cells against how sharply two libraries can
/// be told apart. It trades that only down to the point where they can be told apart at
/// all — pooled outright, one library's fitted rate moves 23 to 38 rungs of the
/// error-rate ladder on nothing but where the search starts (research note §2.2), which
/// is a qualitative break at the bottom of the range rather than a knob position.
///
/// **Four is carried over from the design rather than picked by the measurement**,
/// which finds no reason to prefer it to two: "the attribution bound is a precision knob
/// that is not currently buying precision", and at three reads a bound of two is equally
/// unbiased on 28% fewer cells, neither losing measurable precision against scoring
/// every read against its own library (research note §2.5,
/// `spec/parameter_prepass_generic.md` §1, §12.8). So this is the cheapest saving of
/// attributed cells available if one is ever wanted.
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
/// from them anyway.
///
/// **At one read group an attributed key carries nothing a pooled one does not — but
/// this type does not collapse the arm, because it cannot see how many libraries the
/// sample has.** Entering a single-library sample through
/// [`DepthAltHistogram::add_site`] is the accumulator's job (Milestone C), and it is
/// what keeps the 1,550 of the 1,707 samples in the tomato archive survey that carry
/// one library keyed exactly as they would be with no multi-library machinery at all.
/// It is also what `spec/parameter_prepass_generic.md` §12.6 rests on — the read-group
/// histogram equals the windowed one folded over its windows, cell for cell, on a
/// single-library sample — which would fail if such a sample's sites carried attributed
/// keys in one table and pooled keys in the other.
///
/// **The attributed arm keeps a depth bin and not an exact depth, and that is measured
/// rather than assumed.** Under the adopted ladder the fit's asymptotic bias stays at
/// 0.054 rungs of the error-rate ladder and at worst 0.3% in a genotype frequency —
/// 0.23% in heterozygosity and 0.30% in the homozygous-non-reference rate — against
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
///   it, while `cells()`'s count and the region-sharding associativity of
///   `spec/parameter_prepass_generic.md` §4 — two shards merging to the histogram of
///   their union — would both move.
/// - **Which arm a site takes follows from its counts**, so it is not a choice a call
///   site gets to make differently from its neighbour.
///
/// **The field order below is the sort order.** `Ord` is derived, so a
/// `BTreeMap<SiteKey, _>` comes out bin-major, then by alternative count, then by
/// attribution — which is the order [`DepthAltHistogram::cells`] documents. Reordering
/// the fields to save padding would change a documented public ordering.
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
    ///
    /// **Four inline, matching [`MAX_ATTRIBUTED_ALT_READS`]**, because four read groups
    /// can contribute one alternative read each and that is the widest listing this arm
    /// admits. At two inline — the architecture's illustrative figure — a site spread
    /// over three or four libraries heap-allocates *once per position*, measured at
    /// 1,000 allocations per 1,000 steady-state sites, on the path whose contract is
    /// that it never allocates after a key's first site. The sixteen extra bytes land
    /// only on keys in the sparse map.
    alt_by_group: SmallVec<[(ReadGroupId, u8); 4]>,
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
    /// this site, in any order; the key it produces is canonical. The counts narrow to
    /// one byte on the way in, which the bound above is what makes safe. Three things
    /// pool instead, and each is the same information as an attributed key rather than
    /// less of it:
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
    /// read sets, and summing the duplicates would hide it. **Checked on the listing as
    /// it arrives, before the arm is chosen**, because the arm that pools is where a
    /// doubled count does the damage: a site whose duplicated entries sum above the
    /// bound would otherwise pool quietly, and [`DepthAltHistogram::add_attributed_site`]'s
    /// own cross-check cannot catch it either, since a caller that double-counted a read
    /// set builds both of its arguments from the same reads and the two agree.
    ///
    /// Also if the listing's counts sum past `u32`.
    #[must_use]
    pub fn attributing(depth_bin: DepthBin, alt_by_group: &[(ReadGroupId, u32)]) -> Self {
        // A pairwise scan rather than a sort: the listing holds one entry per read
        // group that showed an alternative read, so it is a handful of entries on a
        // path walked once per covered position, and sorting a copy of it here would
        // put an allocation on that path to save comparisons that cost nothing.
        for (position, &(group, _)) in alt_by_group.iter().enumerate() {
            assert!(
                !alt_by_group[position + 1..]
                    .iter()
                    .any(|&(other, _)| other == group),
                "a site listed the same read group twice: {alt_by_group:?}"
            );
        }

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

        let mut attributed: SmallVec<[(ReadGroupId, u8); 4]> = alt_by_group
            .iter()
            .filter(|&&(_, alt)| alt > 0)
            // Checked at the point of use rather than argued from the early return
            // above and the const assertion at the top of the file. Both hold today; an
            // edit that moves either would otherwise leave a silent truncation here —
            // 256 alternative reads becoming zero — which is neither a compile error
            // nor a panic.
            .map(|&(group, alt)| {
                let attributed = u8::try_from(alt).unwrap_or_else(|_| {
                    panic!(
                        "read group {} showed {alt} alternative reads at a site whose \
                         total is at most {MAX_ATTRIBUTED_ALT_READS}",
                        group.get()
                    )
                });
                (group, attributed)
            })
            .collect();
        attributed.sort_unstable_by_key(|&(group, _)| group);

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

    /// Which library each alternative read came from, or that this key pooled them.
    #[inline]
    #[must_use]
    pub fn attribution(&self) -> Attribution<'_> {
        if self.alt_by_group.is_empty() {
            Attribution::Pooled
        } else {
            Attribution::ByReadGroup(&self.alt_by_group)
        }
    }
}

/// Which library each of a site's alternative reads came from — or that its key does
/// not say.
///
/// **A two-variant view rather than an `Option<&[…]>`, so that a consumer cannot reach
/// for `unwrap_or(&[])`.** An empty listing and a pooled key mean opposite things: the
/// first would say no library showed an alternative read, where the second says several
/// may have and this key forgot which. Collapsing the two gives a wrong likelihood term
/// with no symptom, in the module whose whole thesis is that its wrong numbers have
/// none. The score for a pooled cell sums over the splits the key forgot rather than
/// inventing one (`arch/parameter_prepass_generic.md` §5.1), so the two arms are
/// different arithmetic and a `match` is what forces a consumer to write both.
///
/// It is a public enum where [`SiteKey`] deliberately is not, and the objection that
/// settled that does not transfer: this is *returned*, never constructed by a caller,
/// so a public variant cannot be used to build a key that breaks an invariant.
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub enum Attribution<'a> {
    /// This key does not say which library the alternative reads came from — which is
    /// not the same as saying none did.
    Pooled,
    /// Ascending by read group, carrying no zero entries, and never empty.
    ByReadGroup(&'a [(ReadGroupId, u8)]),
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
/// no symptom. A panic naming the cell is the alternative, and it was measured at
/// 0.2 ns a site — 1.31 ns against a wrapping add's 1.09 over 200 million calls, which
/// is 0.16 s across a whole human genome.
///
/// `Send + Sync` because a shard's tables cross a thread boundary on the way to the
/// merge (`arch/parameter_prepass_generic.md` §3). The seal already guarantees both;
/// stating them here is what saves every helper generic over `C` from repeating the
/// bound, and a bound repeated at many sites is one forgotten at one of them.
pub trait CellCounter:
    Copy + Into<u64> + From<u32> + std::fmt::Debug + Send + Sync + sealed::Sealed
{
    /// Add, or `None` if the sum would not fit. The caller supplies the message,
    /// because what a reader needs is which cell overflowed and not which integer type.
    #[must_use]
    fn checked_add(self, other: Self) -> Option<Self>;
}

impl CellCounter for u32 {
    #[inline]
    fn checked_add(self, other: Self) -> Option<Self> {
        u32::checked_add(self, other)
    }
}

impl CellCounter for u64 {
    #[inline]
    fn checked_add(self, other: Self) -> Option<Self> {
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

/// A tally of what a sample's sites looked like: for each depth bin, one counter per
/// alternative-read count from none to that bin's deepest depth.
///
/// **Dense and ragged.** Ragged because a site's alternative count cannot exceed its
/// depth, so a rectangular table would waste half of itself; dense because iteration
/// order is then fixed, which every fit needs — a fit is a floating-point sum over
/// cells, and floating-point addition is not associative, so a table that walked its
/// cells in a different order between runs would let a fitted rate wobble on identical
/// data. The sites whose alternative reads keep their library sit in a sparse map
/// beside it, because a multi-library sample occupies few of those keys and a
/// single-library one occupies none at all.
///
/// **Holds the binning rule by [`Arc`], because it needs it to address its own
/// cells** — the row offsets live there. A sample has about 8,000 of these, one per
/// 100 kb window, so that is 8,000 pointers to one shared object: 64 kB against the
/// tables' 37 MB. Sharing rather than copying is what lets the fold prove two tables
/// are binned the same way by pointer identity rather than by comparing lengths and
/// hoping (`arch/parameter_prepass_generic.md` §2.2).
///
/// **Generic over its counter width, and the two widths are not interchangeable.** A
/// window accumulates in `u32`; the fold over a sample's windows produces `u64`. The
/// site count alone would not force that: folded over a human genome it reaches
/// 3.1 × 10⁹ against a `u32` ceiling of 4.29 × 10⁹ — close, but inside. **The per-cell
/// depth sum is what forces it**, and it is not close: 3.1 × 10¹¹, seventy-two times
/// the ceiling. A fold that widened the counts and left the depth sums alone would wrap
/// the one quantity `mean_depth_in_cell` exists to hold, which the design records as a
/// mistake that has already shipped once. Those sums arrive with Milestone B3 and the
/// widening with B4; the trait is here now so that neither can be added at the wrong
/// width.
///
/// **The width defaults to `u32`, the window's.** A field or signature written
/// `DepthAltHistogram` with no parameter is a *window* table; only the whole-sample
/// fold spells `DepthAltHistogram<u64>`. Inference will not pick the default up —
/// `DepthAltHistogram::new` with neither a turbofish nor an annotation does not
/// compile — so the default is only ever reached by naming the type, which is where a
/// reader will not see it unless it is said here.
#[derive(Debug)]
pub struct DepthAltHistogram<C: CellCounter = u32> {
    /// The pooled arm: one counter per cell, rows located through
    /// [`DepthBinEdges::row_start`].
    counts: Vec<C>,
    /// The attributed arm: the sites whose alternative reads kept the library each
    /// came from. Sparse, and **empty for a single-library sample provided the
    /// accumulator entered it through [`DepthAltHistogram::add_site`]** (Milestone C) —
    /// this type will populate the map for a lone read group if it is asked to, because
    /// a key cannot see how many libraries a sample has.
    attributed_cells: BTreeMap<SiteKey, C>,
    edges: Arc<DepthBinEdges>,
    /// How many reference positions the loci entered here covered.
    ///
    /// **`u64` rather than `C`, unlike the cells.** The memory argument that forces a
    /// width choice on the cell vectors — four bytes for a site count and, once
    /// Milestone B3 adds the depth sums, four more: eight bytes a cell, ~4.7 kB a
    /// window, ~37 MB a tomato sample (`spec/parameter_prepass_generic.md` §9) — does
    /// not reach one scalar per table, where the whole difference is eight bytes
    /// against four. And the quantity is uncomfortably placed for 32 bits: a read-group
    /// table accumulates over the whole genome without a windowed fold to widen it, and
    /// a human genome's 3.1 × 10⁹ analysable positions sit at 72% of the `u32` ceiling.
    covered_positions: u64,
}

impl<C: CellCounter> DepthAltHistogram<C> {
    /// An empty table binned by `edges`.
    ///
    /// Every table a run creates must be handed the **same** edges object, or their
    /// cells mean different things and folding them together is meaningless.
    #[must_use]
    pub fn new(edges: Arc<DepthBinEdges>) -> Self {
        let counts = vec![C::from(0); edges.cell_count()];
        Self {
            counts,
            attributed_cells: BTreeMap::new(),
            edges,
            covered_positions: 0,
        }
    }

    /// Add one site whose alternative reads are not attributed to a library: every
    /// entry of the read-group histogram, where the whole entry belongs to one group
    /// already, and every entry of a single-library sample.
    ///
    /// The depth bin is derived **here**, from the exact depth the pair carries, and
    /// the depth recorded is the one the site was entered at — the capped one where
    /// the subsampling of Milestone C2 fired. Recording the pre-cap depth would put a
    /// depth in the cell that no site in it actually had.
    ///
    /// `covered` is how many reference positions the locus spanned, which is one for
    /// the overwhelming majority and more where a locus was widened to an indel's
    /// reference span. It is a second argument because a locus's span is not derivable
    /// from its reads, and [`DepthAltHistogram::total_covered_positions`] — what the
    /// inbreeding fit weights windows by — cannot be accumulated without it.
    ///
    /// # Panics
    ///
    /// If the site arrived deeper than the ladder's cap — see
    /// [`DepthAltHistogram::check_within_the_cap`]. If a counter or the covered-position
    /// total runs past its width, which is reported rather than wrapped.
    pub fn add_site(&mut self, site: DepthAndAltReads, covered: Bp) {
        self.check_within_the_cap(site.depth());
        let depth_bin = self.edges.bin_for(site.depth());
        self.add_pooled(depth_bin, site.alt_reads());
        self.add_covered_positions(covered);
    }

    /// Add one site whole, keeping which library each of its alternative reads came
    /// from where that is worth keeping — the windowed histogram of a sample with more
    /// than one read group.
    ///
    /// `alt_by_group` is one entry per read group that showed an alternative read, in
    /// any order. [`SiteKey::attributing`] decides the arm: a site with nothing to
    /// attribute, or with more than [`MAX_ATTRIBUTED_ALT_READS`], joins the dense
    /// table beside the pooled sites rather than the sparse map.
    ///
    /// **Whether a sample's sites come through here or through
    /// [`DepthAltHistogram::add_site`] is a decision about the *sample*, not about the
    /// site.** A sample with one library has nothing to attribute at any site, and
    /// feeding it through here would build a sparse map that
    /// `spec/parameter_prepass_generic.md` §9 does not price and break §12.6's cell-for-cell
    /// equality between the two histograms. That decision belongs to the accumulator
    /// (Milestone C).
    ///
    /// # Panics
    ///
    /// If the breakdown does not account for exactly the site's alternative reads —
    /// the two are counts of the same thing from the same reads, so a disagreement is
    /// one of them being counted from the wrong set. And for every reason
    /// [`SiteKey::attributing`] and [`DepthAltHistogram::add_site`] give, since this
    /// call goes through each: a read group listed twice, an alternative total past
    /// `u32`, a depth above the cap, and a counter past its width.
    pub fn add_attributed_site(
        &mut self,
        site: DepthAndAltReads,
        alt_by_group: &[(ReadGroupId, u32)],
        covered: Bp,
    ) {
        self.check_within_the_cap(site.depth());
        let depth_bin = self.edges.bin_for(site.depth());
        let key = SiteKey::attributing(depth_bin, alt_by_group);
        assert_eq!(
            key.alt_reads(),
            site.alt_reads(),
            "a site showing {} alternative reads was broken down into {} — these are \
             different counts of the same thing, so one of them was counted from the \
             wrong read set",
            site.alt_reads(),
            key.alt_reads()
        );

        match key.attribution() {
            Attribution::ByReadGroup(_) => {
                let counter = self.attributed_cells.entry(key).or_insert(C::from(0));
                Self::count_one_more(counter, depth_bin, site.alt_reads());
            }
            Attribution::Pooled => self.add_pooled(depth_bin, site.alt_reads()),
        }
        self.add_covered_positions(covered);
    }

    /// Every cell that holds at least one site, materialised once.
    ///
    /// **A `Vec` rather than an iterator, because the profile scan re-walks these 161
    /// times** — once per rung of the error-rate ladder — and an iterator would
    /// re-derive the attributed arm's keys on every pass. A few hundred to a few
    /// thousand cells is kilobytes (`arch/parameter_prepass_generic.md` §4.2).
    ///
    /// The ploidy travels with each cell because one error rate is fitted per read
    /// group **across** the ploidies that group covered, so a single scan sees cells of
    /// more than one.
    ///
    /// The order is fixed: the pooled cells bin by bin and, within a bin, by
    /// alternative count, then the attributed ones in key order. Empty cells are left
    /// out — they contribute nothing to any sum, and at three reads a site the ladder
    /// is mostly empty.
    #[must_use]
    pub fn cells(&self, ploidy: Ploidy) -> Vec<Cell> {
        // Both arms' upper bound: at most one cell per ladder cell, plus the sparse
        // map. Over-reserving a few tens of kilobytes once beats the nine reallocations
        // a full table costs otherwise, on a vector the profile scan re-walks 161 times.
        let mut cells = Vec::with_capacity(self.counts.len() + self.attributed_cells.len());

        for depth_bin in self.edges.bins() {
            let row = self.edges.row_start(depth_bin);
            let width = *self.edges.depth_range(depth_bin).end() as usize + 1;
            for (alt_reads, count) in self.counts[row..row + width].iter().enumerate() {
                let sites: u64 = (*count).into();
                if sites > 0 {
                    cells.push(Cell {
                        key: SiteKey::pooled(depth_bin, alt_reads as u32),
                        ploidy,
                        sites,
                    });
                }
            }
        }
        cells.extend(self.attributed_cells.iter().map(|(key, count)| Cell {
            key: key.clone(),
            ploidy,
            sites: (*count).into(),
        }));
        cells
    }

    /// How many **loci** entered — not how many reference positions they covered.
    ///
    /// **Derived, never stored** (`spec/parameter_prepass_generic.md` §4): every locus
    /// enters exactly one cell, including the overwhelming majority that show no
    /// alternative read, so the cell counts already carry it and a second counter
    /// could only disagree with them.
    ///
    /// # Panics
    ///
    /// If the cells sum past `u64`, which they cannot through this module's own entry
    /// points — 583 `u32` cells cap the sum at about 2.5 × 10¹² — so the check is
    /// against a later fold, not against a reachable input.
    #[must_use]
    pub fn total_loci(&self) -> u64 {
        let cells = self.counts.iter().chain(self.attributed_cells.values());
        cells.fold(0u64, |total, &count| {
            let in_cell: u64 = count.into();
            total.checked_add(in_cell).unwrap_or_else(|| {
                panic!(
                    "this table's locus count passed u64 at {total}, with a cell holding \
                     {in_cell} more"
                )
            })
        })
    }

    /// How many reference **positions** those loci covered — `Σ region.len()`.
    ///
    /// Not derivable from [`DepthAltHistogram::total_loci`], because a generic locus
    /// widened to an indel's reference span is one entry and several positions. This
    /// is what the runs model weights window posteriors by, so that the inbreeding
    /// coefficient is a fraction of the analysable genome rather than of the locus
    /// list and a window dense in indels is not under-weighted.
    #[must_use]
    pub fn total_covered_positions(&self) -> u64 {
        self.covered_positions
    }

    /// Refuse a site the ladder has no bin for.
    ///
    /// **[`DepthBinEdges::bin_for`] is total by design and answers the top bin for any
    /// depth**, which is right for a binning rule and wrong as a place to *enter* such a
    /// site: a 500-read site recorded in the 98–124 bin puts a depth in that cell that
    /// no site in it could have had, and Milestone B3's per-cell depth sum is exactly
    /// the quantity that would carry it. That is the same class of wrong number
    /// `mean_depth_in_cell` exists to prevent, arriving from above instead of below.
    ///
    /// The cap is enforced upstream, by the subsampling of Milestone C2 — this makes
    /// that a checked precondition rather than an assumed one, for one comparison on a
    /// path that already looks the depth up.
    fn check_within_the_cap(&self, depth: u32) {
        let cap = self.edges.max_depth();
        assert!(
            depth <= cap,
            "a site of depth {depth} reached a histogram whose ladder tops out at \
             {cap} — a deeper site has to be subsampled down to the cap before it is \
             entered, not binned into the top bin"
        );
    }

    /// Where a pooled cell sits in the flat vector.
    ///
    /// **The bound check is the one that matters in this file.** A bin's row is only as
    /// wide as its deepest depth, so an alternative count above that would address the
    /// *next* bin's cells — a real counter, incremented, in a table that still sums to
    /// the right total. Nothing downstream would show it.
    ///
    /// It is unreachable through the two entry points, and deliberately kept: with
    /// `alt_reads <= depth` from [`DepthAndAltReads`] and `depth <= cap` from
    /// [`DepthAltHistogram::check_within_the_cap`], a site's alternative count cannot
    /// exceed its own bin's top. What it guards is the *other* way in — Milestone B3
    /// looks a cell up by a [`SiteKey`] its caller supplies, and a key naming a bin and
    /// an alternative count that do not go together has no such guarantee behind it.
    fn pooled_cell_index(&self, depth_bin: DepthBin, alt_reads: u32) -> usize {
        let deepest = *self.edges.depth_range(depth_bin).end();
        assert!(
            alt_reads <= deepest,
            "a site in depth bin {} showed {alt_reads} alternative reads, above the \
             bin's deepest depth {deepest} — a site deeper than the ladder's cap has to \
             be subsampled down to it before it reaches a histogram",
            depth_bin.get()
        );
        self.edges.row_start(depth_bin) + alt_reads as usize
    }

    /// Count one more site into the dense table. The two entry points share it so that
    /// the bound check and the increment cannot come apart, and so that Milestone B4's
    /// fold has one place to add through.
    fn add_pooled(&mut self, depth_bin: DepthBin, alt_reads: u32) {
        let index = self.pooled_cell_index(depth_bin, alt_reads);
        Self::count_one_more(&mut self.counts[index], depth_bin, alt_reads);
    }

    fn count_one_more(counter: &mut C, depth_bin: DepthBin, alt_reads: u32) {
        *counter = counter.checked_add(C::from(1)).unwrap_or_else(|| {
            panic!(
                "the cell for depth bin {} at {alt_reads} alternative reads holds more \
                 sites than a {}-bit counter can, so this table needs the whole-sample \
                 fold's width",
                depth_bin.get(),
                std::mem::size_of::<C>() * 8
            )
        });
    }

    fn add_covered_positions(&mut self, covered: Bp) {
        let so_far = self.covered_positions;
        self.covered_positions = so_far.checked_add(covered.get()).unwrap_or_else(|| {
            panic!(
                "this table's covered positions passed u64 at {so_far}, with a locus \
                 spanning {} more",
                covered.get()
            )
        });
    }
}

/// One cell of the tally, as a fit reads it.
///
/// **A struct rather than the `(SiteKey, Ploidy, u64)` the architecture sketches**
/// (`arch/parameter_prepass_generic.md` §2.2), because the third member is a count of
/// sites and nothing in a tuple says which of the three that is. This vector is the seam
/// between `generic/` and `fitting/` — the one interface the project's shaping-versus-
/// mathematics split exists to keep clean — and the profile scan walks it 161 times per
/// fit, in code that would otherwise be written in `.2`. Introduced now because
/// retro-fitting it once the fits exist would touch every one of them.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Cell {
    /// What the sites in this cell looked like.
    pub key: SiteKey,
    /// The genotype set this cell is to be scored against. It travels with the cell
    /// because one error rate is fitted per read group **across** the ploidies that
    /// group covered, so a single scan sees cells of more than one.
    pub ploidy: Ploidy,
    /// How many sites landed here.
    pub sites: u64,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn group(id: u32) -> ReadGroupId {
        ReadGroupId(id)
    }

    fn ladder() -> Arc<DepthBinEdges> {
        Arc::new(DepthBinEdges::new())
    }

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a positive copy number")
    }

    /// One reference position, which is what all but the widened loci cover.
    const ONE_POSITION: Bp = Bp(1);

    /// A depth and an alternative count are two `u32`s in the same order of magnitude,
    /// so the pair is the one place that can tell them apart — and it does, by
    /// refusing the transposition rather than storing it.
    #[test]
    fn a_depth_and_an_alternative_count_keep_the_order_they_were_given_in() {
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
            Attribution::ByReadGroup(&[(group(0), 1u8), (group(7), 2u8)]),
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
        assert_eq!(
            listed.attribution(),
            Attribution::ByReadGroup(&[(group(0), 2u8)])
        );
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
        assert_eq!(nothing.attribution(), Attribution::Pooled);

        let at_the_bound = SiteKey::attributing(bin, &[(group(0), 2), (group(7), 2)]);
        assert_eq!(at_the_bound.alt_reads(), MAX_ATTRIBUTED_ALT_READS);
        assert!(
            matches!(at_the_bound.attribution(), Attribution::ByReadGroup(_)),
            "the bound is inclusive, so four alternative reads keep their libraries"
        );

        let over_the_bound = SiteKey::attributing(bin, &[(group(0), 2), (group(7), 3)]);
        assert_eq!(over_the_bound, SiteKey::pooled(bin, 5));
        assert_eq!(over_the_bound.attribution(), Attribution::Pooled);
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

    /// The key sorts **bin-major**, then by alternative count, then by attribution —
    /// which is what makes a `BTreeMap` of these iterate in an order that does not vary
    /// between runs, and every fit is a floating-point sum over cells, which is not
    /// associative.
    ///
    /// **The bin-7-alt-9 key is what makes this test able to fail.** `Ord` is derived,
    /// so the order is carried by the order the fields are declared in, and with every
    /// key's bin and alternative count ascending together an alt-major sort produces the
    /// identical vector — swapping the two fields left the whole module green. This one
    /// key has them disagree: bin 9 with 7 alternative reads must come *before* bin 12
    /// with 3.
    #[test]
    fn keys_order_by_depth_bin_before_alternative_count() {
        let mut keys = [
            SiteKey::pooled(DepthBin(12), 3),
            SiteKey::pooled(DepthBin(9), 7),
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
                    matches!(key.attribution(), Attribution::ByReadGroup(_)),
                )
            })
            .collect();
        assert_eq!(
            shape,
            vec![
                (9, 0, false),
                (9, 1, false),
                (9, 1, true),
                (9, 7, false),
                (12, 3, false),
            ]
        );
    }

    /// The bound is stated as a literal rather than left to be read off a cell count,
    /// because it is a number a later reader will want to change — and the research
    /// note's answer is that two would do as well, not that four was chosen.
    #[test]
    fn the_attribution_bound_is_four() {
        assert_eq!(MAX_ATTRIBUTED_ALT_READS, 4);
    }

    /// **A counter that ran past its width must panic, not wrap.** The release profile
    /// this repo builds leaves `overflow-checks` off, so a wrapped depth sum would come
    /// back as a small number and be scored as a shallow cell — which is the failure
    /// the per-cell depth sum exists to prevent, arriving by another route.
    #[test]
    fn a_counter_reports_the_sum_it_cannot_hold_rather_than_wrapping() {
        assert_eq!(CellCounter::checked_add(1u32, 2u32), Some(3u32));
        assert_eq!(CellCounter::checked_add(u32::MAX, 1u32), None);

        // The width the whole-sample fold widens to holds what `u32` cannot: a human
        // genome's depth sum reaches 3.1 × 10¹¹ against a `u32` ceiling of 4.29 × 10⁹.
        let human_depth_sum = 310_000_000_000u64;
        assert!(human_depth_sum > u64::from(u32::MAX));
        assert_eq!(
            CellCounter::checked_add(human_depth_sum, 1),
            Some(310_000_000_001)
        );
        assert_eq!(CellCounter::checked_add(u64::MAX, 1u64), None);
    }

    /// Both widths start empty and widen upward without loss, which is what the fold
    /// of Milestone B4 rests on. An empty cell is spelled `C::from(0)` through the
    /// trait's own `From<u32>` supertrait rather than a `ZERO` associated constant,
    /// which would be six lines for a spelling.
    #[test]
    fn both_counter_widths_start_at_zero_and_widen_to_u64() {
        fn empty_cell<C: CellCounter>() -> u64 {
            C::from(0).into()
        }

        assert_eq!(empty_cell::<u32>(), 0);
        assert_eq!(empty_cell::<u64>(), 0);
        assert_eq!(Into::<u64>::into(u32::MAX), 4_294_967_295u64);
        assert_eq!(<u64 as From<u32>>::from(7), 7u64);
    }

    /// **Every cell of the ladder is reachable, and no two sites share one — which is
    /// what says the rows are the right widths.** A site is entered at each bin's
    /// deepest depth once per legal alternative count, so the 583 writes should land on
    /// 583 distinct cells holding one site each. A row one cell too narrow or a row
    /// offset one cell out puts two of those writes in the same counter and leaves
    /// another at zero, which this reads off directly: the cell count drops and a cell
    /// comes back holding two.
    ///
    /// It is the only test here that would notice a mis-sized row, because every other
    /// one uses a handful of cells that a wrong offset could still keep apart.
    #[test]
    fn every_cell_of_the_ladder_holds_exactly_the_sites_addressed_to_it() {
        let edges = ladder();
        let mut table = DepthAltHistogram::<u32>::new(edges.clone());
        let mut written = 0u64;

        for bin in 0..edges.bin_count() {
            let depth_bin = DepthBin(bin as u16);
            let deepest = *edges.depth_range(depth_bin).end();
            for alt_reads in 0..=deepest {
                table.add_site(DepthAndAltReads::new(deepest, alt_reads), ONE_POSITION);
                written += 1;
            }
        }

        assert_eq!(
            written, 583,
            "the ladder's own cell count (research note §4.3)"
        );
        let cells = table.cells(diploid());
        assert_eq!(
            cells.len(),
            583,
            "each (bin, alternative count) pair is its own cell"
        );
        assert!(
            cells.iter().all(|cell| cell.sites == 1),
            "two sites landed in one cell, so a row is mis-sized or mis-placed"
        );
        assert_eq!(table.total_loci(), written);
    }

    /// **A site deeper than the ladder's cap is refused on both arms**, rather than
    /// binned into the top bin. `bin_for` is total by design and answers bin 19 for any
    /// depth, so without this the entry points accepted a 500-read site silently — and
    /// Milestone B3's per-cell depth sum would then hold a depth no site in that bin
    /// could have had.
    ///
    /// The two arms are tested separately because the attributed one does not go
    /// through the dense table's addressing at all, so it had no depth check of any
    /// kind: a 500-read site with two alternative reads went straight into the sparse
    /// map.
    #[test]
    #[should_panic(expected = "ladder tops out at 124")]
    fn a_site_deeper_than_the_cap_is_refused_on_the_pooled_arm() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_site(DepthAndAltReads::new(500, 5), ONE_POSITION);
    }

    #[test]
    #[should_panic(expected = "ladder tops out at 124")]
    fn a_site_deeper_than_the_cap_is_refused_on_the_attributed_arm() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_attributed_site(
            DepthAndAltReads::new(500, 2),
            &[(group(0), 1), (group(1), 1)],
            ONE_POSITION,
        );
    }

    /// The row bound, reached the only way it now can be: by naming a bin and an
    /// alternative count that do not go together. The two entry points cannot — a
    /// site's alternative count is at most its depth and its depth is at most the cap —
    /// but Milestone B3 looks a cell up by a key its caller supplies, and an
    /// alternative count above the bin's top would address the *next* bin's cells: a
    /// real counter, in a table that still sums to the right total.
    #[test]
    #[should_panic(expected = "above the bin's deepest depth 5")]
    fn a_cell_named_outside_its_bins_row_is_refused() {
        let table = DepthAltHistogram::<u32>::new(ladder());
        let _ = table.pooled_cell_index(DepthBin(5), 99);
    }

    /// **The cells come back in one order however the sites arrived**, which every fit
    /// depends on: a fit is a floating-point sum over cells, and floating-point
    /// addition is not associative, so an order that varied would let a fitted rate
    /// wobble between runs of identical data.
    #[test]
    fn the_cells_come_back_in_one_order_whatever_order_the_sites_arrived_in() {
        let from_seven = [(group(7), 1u32)];
        let from_zero = [(group(0), 2u32)];
        let sites: [(DepthAndAltReads, &[(ReadGroupId, u32)]); 5] = [
            (DepthAndAltReads::new(30, 15), &[]),
            (DepthAndAltReads::new(12, 1), &from_seven),
            (DepthAndAltReads::new(3, 0), &[]),
            (DepthAndAltReads::new(12, 2), &from_zero),
            (DepthAndAltReads::new(30, 1), &[]),
        ];

        fn fill<'a>(
            order: impl Iterator<Item = &'a (DepthAndAltReads, &'a [(ReadGroupId, u32)])>,
        ) -> Vec<Cell> {
            let mut table = DepthAltHistogram::<u32>::new(ladder());
            for &(site, by_group) in order {
                if by_group.is_empty() {
                    table.add_site(site, ONE_POSITION);
                } else {
                    table.add_attributed_site(site, by_group, ONE_POSITION);
                }
            }
            table.cells(diploid())
        }

        let forwards = fill(sites.iter());
        let backwards = fill(sites.iter().rev());
        assert_eq!(forwards, backwards);
        assert_eq!(forwards.len(), 5);

        // Bin-major, then by alternative count, then the attributed cells. Depths 11
        // to 13 share bin 10 and 29 to 36 share bin 14, so the two depth-12 sites and
        // the two depth-30 ones are neighbours here rather than far apart.
        let shape: Vec<(u16, u32, bool)> = forwards
            .iter()
            .map(|cell| {
                (
                    cell.key.depth_bin().get(),
                    cell.key.alt_reads(),
                    matches!(cell.key.attribution(), Attribution::ByReadGroup(_)),
                )
            })
            .collect();
        assert_eq!(
            shape,
            vec![
                (3, 0, false),
                (14, 1, false),
                (14, 15, false),
                (10, 1, true),
                (10, 2, true),
            ]
        );
    }

    /// **The two counters are not the same number, and the difference is the whole
    /// reason the second one is accumulated.** A generic locus widened to an indel's
    /// reference span is one entry in the tally and several covered positions, and it
    /// is the covered positions the inbreeding fit weights windows by.
    #[test]
    fn a_locus_widened_to_an_indels_span_is_one_entry_and_several_positions() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());

        table.add_site(DepthAndAltReads::new(12, 0), ONE_POSITION);
        table.add_site(DepthAndAltReads::new(12, 1), Bp(4));

        assert_eq!(table.total_loci(), 2);
        assert_eq!(table.total_covered_positions(), 5);
    }

    /// **The same accounting on the attributed arm, which is the entry point for every
    /// window of a multi-library sample.** Deleting that arm's span accumulation left
    /// the whole suite green: every other attributed test asserts on `cells()` or
    /// `total_loci()`, and neither reads the covered positions. Lost, every such window
    /// would report zero covered positions, the runs model would weight them all at
    /// zero, and the inbreeding coefficient would come out of a weighting that nothing
    /// panics on and nothing prints.
    #[test]
    fn a_widened_locus_counts_its_whole_span_on_the_attributed_arm_too() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());

        table.add_site(DepthAndAltReads::new(12, 0), ONE_POSITION);
        table.add_attributed_site(DepthAndAltReads::new(12, 1), &[(group(0), 1)], Bp(4));

        assert_eq!(table.total_loci(), 2);
        assert_eq!(
            table.total_covered_positions(),
            5,
            "an attributed locus's span counts exactly as a pooled one's does"
        );
    }

    /// A span total past `u64` is reported rather than wrapped. [`Bp`] is a public
    /// newtype over `u64`, so a caller can hand this any span, and a wrapped total would
    /// report a near-zero covered-position count for a whole sample — the weight the
    /// inbreeding fit divides by.
    #[test]
    #[should_panic(expected = "covered positions passed u64")]
    fn a_covered_position_total_past_u64_is_reported_rather_than_wrapped() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_site(DepthAndAltReads::new(1, 0), Bp(u64::MAX));
        table.add_site(DepthAndAltReads::new(1, 0), Bp(1));
    }

    /// **The wiring, not just the trait.** `CellCounter::checked_add` is pinned above,
    /// but what `count_one_more` does with a `None` was not: swapping its panic for
    /// `unwrap_or(*counter)` — a silently dropped increment, which is the outcome the
    /// trait exists to prevent — left the suite green. Reached directly, because
    /// filling a `u32` cell through `add_site` would take 4.29 × 10⁹ calls.
    #[test]
    #[should_panic(expected = "holds more sites than a 32-bit counter can")]
    fn a_cell_at_its_counters_ceiling_names_itself_rather_than_wrapping() {
        let mut counter = u32::MAX;
        DepthAltHistogram::<u32>::count_one_more(&mut counter, DepthBin(3), 2);
    }

    /// The width the whole-sample fold produces counts and reports exactly as the
    /// window width does, which is what Milestone B4's widening rests on — and until
    /// this test, `DepthAltHistogram<u64>` was never instantiated at all.
    #[test]
    fn the_wide_table_counts_and_reports_exactly_as_the_narrow_one_does() {
        let mut wide = DepthAltHistogram::<u64>::new(ladder());
        let mut narrow = DepthAltHistogram::<u32>::new(ladder());

        for table in [&mut narrow as &mut dyn FilledLikeTheOther, &mut wide] {
            table.fill();
        }

        assert_eq!(wide.cells(diploid()), narrow.cells(diploid()));
        assert_eq!(wide.total_loci(), narrow.total_loci());
        assert_eq!(
            wide.total_covered_positions(),
            narrow.total_covered_positions()
        );
        assert_eq!(wide.total_loci(), 2);
    }

    /// Lets one filling recipe drive both counter widths, so the two tables cannot be
    /// filled differently by accident and the comparison above stay green.
    trait FilledLikeTheOther {
        fn fill(&mut self);
    }

    impl<C: CellCounter> FilledLikeTheOther for DepthAltHistogram<C> {
        fn fill(&mut self) {
            self.add_site(DepthAndAltReads::new(12, 1), ONE_POSITION);
            self.add_attributed_site(DepthAndAltReads::new(12, 1), &[(group(0), 1)], Bp(4));
        }
    }

    /// The architecture's worked example, run through the ladder that actually bins
    /// it. Ten positions become **three** cells and not the six the illustration
    /// shows, because depths 11 to 13 share one bin and 29 to 36 share another — the
    /// illustration is drawn on exact depths, which is what the ladder replaces.
    #[test]
    fn ten_positions_become_three_cells_once_their_depths_are_binned() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        for site in [
            (12, 0),
            (12, 0),
            (11, 0),
            (12, 1),
            (13, 0),
            (12, 0),
            (30, 15),
            (12, 0),
            (11, 1),
            (12, 0),
        ] {
            table.add_site(DepthAndAltReads::new(site.0, site.1), ONE_POSITION);
        }

        let cells: Vec<(u16, u32, u64)> = table
            .cells(diploid())
            .into_iter()
            .map(|cell| (cell.key.depth_bin().get(), cell.key.alt_reads(), cell.sites))
            .collect();
        assert_eq!(cells, vec![(10, 0, 7), (10, 1, 2), (14, 15, 1)]);
        assert_eq!(table.total_loci(), 10);
    }

    /// Two sites with the same depth and the same alternative count are **different
    /// cells** when one of them says which library its alternative reads came from.
    /// That is the whole point of the attributed arm: without it the likelihood is
    /// exactly flat along every combination of per-library rates holding their
    /// share-weighted mean fixed.
    #[test]
    fn an_attributed_site_and_a_pooled_one_of_the_same_shape_are_separate_cells() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());

        table.add_attributed_site(
            DepthAndAltReads::new(20, 2),
            &[(group(1), 1), (group(0), 1)],
            ONE_POSITION,
        );
        table.add_site(DepthAndAltReads::new(20, 2), ONE_POSITION);

        let cells = table.cells(diploid());
        assert_eq!(cells.len(), 2);
        assert_eq!(
            cells[0].key.attribution(),
            Attribution::Pooled,
            "the dense table comes first"
        );
        assert_eq!(
            cells[1].key.attribution(),
            Attribution::ByReadGroup(&[(group(0), 1u8), (group(1), 1u8)])
        );
        assert_eq!(table.total_loci(), 2);
    }

    /// A multi-library site showing no alternative read has no attribution to keep, so
    /// it joins the dense table — which is where it belongs, since that cell holds the
    /// overwhelming majority of the genome and the sparse map is sized for the few.
    #[test]
    fn a_multi_library_site_with_no_alternative_read_joins_the_dense_table() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_attributed_site(
            DepthAndAltReads::new(20, 0),
            &[(group(0), 0), (group(1), 0)],
            ONE_POSITION,
        );

        let cells = table.cells(diploid());
        assert_eq!(cells.len(), 1);
        assert_eq!(cells[0].key, SiteKey::pooled(DepthBin(12), 0));
        assert_eq!(cells[0].sites, 1);
    }

    #[test]
    #[should_panic(expected = "different counts of the same thing")]
    fn a_breakdown_that_does_not_account_for_the_sites_alternative_reads_is_refused() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_attributed_site(DepthAndAltReads::new(20, 3), &[(group(0), 1)], ONE_POSITION);
    }

    /// An empty table reports emptiness rather than a row of zeros: a cell nothing
    /// landed in contributes nothing to any sum, and at three reads a site most of the
    /// ladder is empty.
    #[test]
    fn an_empty_table_holds_no_cells_and_counts_nothing() {
        let table = DepthAltHistogram::<u32>::new(ladder());

        assert!(table.cells(diploid()).is_empty());
        assert_eq!(table.total_loci(), 0);
        assert_eq!(table.total_covered_positions(), 0);
    }

    /// The ploidy a scan will score each cell against travels with the cell, because
    /// one error rate is fitted per read group across every ploidy that group covered
    /// — so a single scan sees cells of more than one.
    #[test]
    fn every_cell_carries_the_ploidy_it_is_to_be_scored_at() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_site(DepthAndAltReads::new(12, 1), ONE_POSITION);
        table.add_site(DepthAndAltReads::new(30, 2), ONE_POSITION);

        let tetraploid = Ploidy::try_new(4).expect("a positive copy number");
        assert!(
            table
                .cells(tetraploid)
                .iter()
                .all(|cell| cell.ploidy == tetraploid)
        );
    }

    /// **`Pooled` and an empty listing mean opposite things, and the view is what keeps
    /// them apart.** A pooled key with five alternative reads says several libraries may
    /// have contributed and this key forgot which; an empty listing would say none did.
    /// Collapsing the two — which an `Option<&[…]>` invites with `unwrap_or(&[])` —
    /// gives a wrong likelihood term with nothing to show it, and the score for a pooled
    /// cell sums over the splits the key forgot rather than inventing one.
    #[test]
    fn a_pooled_key_is_not_an_empty_listing() {
        let bin = DepthBin(19);
        let pooled = SiteKey::attributing(bin, &[(group(0), 2), (group(7), 3)]);

        assert_eq!(pooled.attribution(), Attribution::Pooled);
        assert_ne!(
            pooled.attribution(),
            Attribution::ByReadGroup(&[]),
            "five alternative reads whose libraries are unknown is not zero of them"
        );
        assert_eq!(pooled.alt_reads(), 5);

        // And the arm that does say: never empty, so `ByReadGroup(&[])` is a state no
        // key produces.
        let attributed = SiteKey::attributing(bin, &[(group(7), 1)]);
        assert_eq!(
            attributed.attribution(),
            Attribution::ByReadGroup(&[(group(7), 1u8)])
        );
    }

    /// **The duplicate guard has to fire on the arm that pools, which is where a
    /// doubled count does the damage.** Before it moved ahead of the arm choice it sat
    /// after the early return, so a listing whose duplicated entries summed above
    /// [`MAX_ATTRIBUTED_ALT_READS`] pooled quietly — and
    /// `add_attributed_site`'s own cross-check cannot catch that either, because a
    /// caller that double-counted a read set builds both of its arguments from the same
    /// reads and the two agree.
    #[test]
    #[should_panic(expected = "listed the same read group twice")]
    fn a_duplicate_listing_that_would_pool_is_refused_too() {
        let _ = SiteKey::attributing(DepthBin(12), &[(group(4), 3), (group(4), 3)]);
    }

    /// The same guard where one of the duplicated entries is a zero, which the
    /// zero-dropping filter used to hide.
    #[test]
    #[should_panic(expected = "listed the same read group twice")]
    fn a_duplicate_listing_with_a_zero_entry_is_refused_too() {
        let _ = SiteKey::attributing(DepthBin(12), &[(group(4), 0), (group(4), 2)]);
    }

    /// A listing whose counts sum past `u32` names the group that pushed it over,
    /// rather than wrapping to a small total and keying the site as a quiet one.
    #[test]
    #[should_panic(expected = "pushed a site's alternative count past u32")]
    fn a_listing_summing_past_u32_names_the_group_that_pushed_it_over() {
        let _ = SiteKey::attributing(DepthBin(12), &[(group(0), u32::MAX), (group(1), 1)]);
    }

    /// **Three groups, in every order, because two cannot tell a sort from a swap.**
    /// The canonical-listing test above uses two entries, where a single "swap if out of
    /// order" pass is indistinguishable from a sort — and wrong on three.
    #[test]
    fn three_read_groups_reach_one_key_from_any_order_they_arrive_in() {
        let bin = DepthBin(12);
        let canonical = SiteKey::attributing(bin, &[(group(0), 1), (group(3), 1), (group(9), 1)]);

        for listing in [
            [(group(9), 1u32), (group(0), 1), (group(3), 1)],
            [(group(3), 1), (group(9), 1), (group(0), 1)],
            [(group(9), 1), (group(3), 1), (group(0), 1)],
            [(group(0), 1), (group(9), 1), (group(3), 1)],
        ] {
            assert_eq!(
                SiteKey::attributing(bin, &listing),
                canonical,
                "{listing:?}"
            );
        }
        assert_eq!(
            canonical.attribution(),
            Attribution::ByReadGroup(&[(group(0), 1u8), (group(3), 1u8), (group(9), 1u8)])
        );
    }

    /// Two attributed cells that tie on depth bin and alternative count are ordered by
    /// their listing — the last field, and the one the other ordering tests never reach
    /// because they never present a tie. Without it, `cells()`'s order across the
    /// attributed arm would depend on insertion history, and a fitted rate would wobble
    /// between runs over identical data.
    #[test]
    fn two_attributed_cells_of_the_same_shape_order_by_their_listing() {
        // Depth 20 falls in bin 12, which is the bin the two keys below are built at.
        let bin = DepthBin(12);
        let from_early = SiteKey::attributing(bin, &[(group(1), 2)]);
        let from_late = SiteKey::attributing(bin, &[(group(8), 2)]);

        assert_eq!(from_early.depth_bin(), from_late.depth_bin());
        assert_eq!(from_early.alt_reads(), from_late.alt_reads());
        assert!(
            from_early < from_late,
            "a tie on bin and count falls through to the read-group listing"
        );

        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_attributed_site(DepthAndAltReads::new(20, 2), &[(group(8), 2)], ONE_POSITION);
        table.add_attributed_site(DepthAndAltReads::new(20, 2), &[(group(1), 2)], ONE_POSITION);
        assert_eq!(
            table
                .cells(diploid())
                .into_iter()
                .map(|cell| cell.key)
                .collect::<Vec<_>>(),
            vec![from_early, from_late],
            "insertion order does not survive into the cells"
        );
    }

    /// The whole path, so the guard is pinned where a caller would actually meet it:
    /// a double-counted read set reaching a cell is a heterozygote invented out of
    /// nothing.
    #[test]
    #[should_panic(expected = "listed the same read group twice")]
    fn a_double_counted_read_set_cannot_reach_a_cell() {
        let mut table = DepthAltHistogram::<u32>::new(ladder());
        table.add_attributed_site(
            DepthAndAltReads::new(20, 6),
            &[(group(4), 3), (group(4), 3)],
            ONE_POSITION,
        );
    }
}
