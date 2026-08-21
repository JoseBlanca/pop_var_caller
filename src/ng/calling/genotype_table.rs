//! Which genotypes a locus is called over, and the three flat tables every step of the
//! calling loop reads them from.
//!
//! A locus's candidate alleles fix its candidate genotypes: at ploidy 2 over alleles
//! `A` and `B` there are three — `AA`, `AB`, `BB` — and at ploidy 4 over six alleles
//! there are 126. Enumerating them, counting each allele's copies in each, and
//! computing how many orderings of the genome's copies spell each one are all pure
//! functions of the pair `(ploidy, allele count)`, and that pair repeats across
//! millions of loci in a run. So the tables are built once per pair and shared
//! ([`GenotypeTable::build`]), which is the port of production's `GenotypeShape`
//! (`src/var_calling/posterior_engine/shape.rs`, `doc/devel/ng/arch/calling_em_loop.md`
//! §8).
//!
//! **The genotype order is the VCF one, and it is a contract rather than an
//! implementation detail.** `PL` and `GL` fields are written in it, so a table that
//! enumerated differently would produce records whose likelihoods name the wrong
//! genotypes — with nothing crashing. [`GenotypeTable::build`] reproduces the order
//! production's `genotype_order` produces; the test module pins it two ways — the
//! diploid triallelic and tetraploid biallelic orders written out by hand, and the
//! ordering rule itself as a law over a grid of shapes — and the value-for-value
//! comparison against production's own function is step C2 of
//! `doc/devel/ng/impl_plan/calling_foundations.md`.

use std::cell::RefCell;
use std::sync::Arc;

use crate::ng::types::{AlleleId, Ploidy};

/// The most alleles one locus can be called over: an [`AlleleId`] names exactly this
/// many, ids 0 through 65,535. A table one allele wider could not say which allele its
/// last row is homozygous for.
pub const MAX_ALLELE_COUNT: usize = u16::MAX as usize + 1;

/// Cache the tables for ploidy up to 8, ported from production's bound
/// (`src/var_calling/posterior_engine/shape.rs`). A diploid or tetraploid locus falls
/// inside; so does the hexaploid or octoploid region of a polyploid crop.
///
/// **Outside the bounds the table is still built, just never kept**: every
/// [`GenotypeTable::build`] for a dodecaploid, or for a locus with 20 candidate
/// alleles, computes a fresh one and hands back a fresh `Arc`. The values are the
/// same; the sharing is what is lost, and past the bound that is not a small loss —
/// at ploidy 8, sixteen alleles is 490,314 genotypes and seventeen is 735,471, so
/// every locus past the bound rebuilds a table of that size instead of copying a
/// pointer.
///
/// **What the bound costs when it is kept.** The slot array itself is 153 pointers,
/// about 1.2 kB a thread. What each filled slot retains is its table, and nothing
/// evicts: the widest cached shape, ploidy 8 over 16 alleles, is 490,314 genotypes
/// and 37,263,864 bytes — 490,314 × (16 × 4 + 8 + 4) — held until the thread exits.
/// A thread that met every shape inside the bounds would hold about 136 MB. That is
/// the deliberate trade, and it is the reason to think twice before widening the
/// bound rather than the slots.
pub const MAX_CACHED_PLOIDY: usize = 8;
/// The allele half of the cache bound; see [`MAX_CACHED_PLOIDY`] for the trade. Sixteen
/// is generous against what a run asks for: the candidate cap ships at 6
/// (`DEFAULT_MAX_CANDIDATE_ALLELES`, `doc/devel/ng/arch/calling_em_loop.md` §8).
pub const MAX_CACHED_ALLELE_COUNT: usize = 16;
/// One slot per `(ploidy, allele count)` pair inside the bounds, both counted from
/// zero so the pair indexes the array directly. Twenty-five of the 153 are unreachable
/// — a [`Ploidy`] cannot be zero and a zero-allele table is refused — and the `+ 1`
/// terms are kept because they make the index arithmetic obvious.
const CACHE_SLOTS: usize = (MAX_CACHED_PLOIDY + 1) * (MAX_CACHED_ALLELE_COUNT + 1);

thread_local! {
    /// Per-thread, so the loop's workers never contend for it and no lock sits on the
    /// path a locus takes. Each thread pays its own first build of a shape; a cohort
    /// run repeats one shape across millions of loci, so that is one build against
    /// millions of hits.
    static TABLE_CACHE: RefCell<[Option<Arc<GenotypeTable>>; CACHE_SLOTS]> =
        const { RefCell::new([const { None }; CACHE_SLOTS]) };
}

/// Which candidate genotype: a row index into one locus's [`GenotypeTable`].
///
/// **The loop's working currency.** A pass scores every sample against every candidate
/// genotype, and carrying a row number rather than an owned multiset of allele ids is
/// what keeps the per-pass work free of allocation
/// (`doc/devel/ng/arch/calling_em_loop.md` §2). The owned [`Genotype`] a sample is
/// finally called as is minted from the winning row on the last pass only.
///
/// **An index means nothing without the table it was minted against**, exactly as an
/// [`AlleleId`] means nothing without its locus: row 4 of a triallelic diploid table is
/// `1/2`, and row 4 of a tetraploid one is something else. The type carries no shape,
/// so an out-of-range index is caught where the table is read — which is why the
/// lookups on [`GenotypeTable`] return `Option` rather than indexing.
///
/// `u32` because the count grows fast with ploidy: 16 alleles at ploidy 8 is 490,314
/// genotypes, already past `u16`, and a `u32` covers every shape that could be built
/// before the table itself exhausted memory.
///
/// [`Genotype`]: crate::ng::types::Genotype
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct GenotypeIdx(pub u32);

impl GenotypeIdx {
    /// The row number itself.
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

/// The candidate genotypes of one `(ploidy, allele count)` shape, with the three
/// quantities the calling loop reads for each of them.
///
/// Built once per shape and shared behind an `Arc` ([`Self::build`]); nothing here
/// depends on a particular locus, only on how many genome copies a sample has and how
/// many alleles the locus is called over.
///
/// The three tables, all in genotype order:
///
/// - **allele counts** — how many copies of each allele the genotype carries. Row *g*
///   of a `genotype_count × allele_count` table, stored flat and row-major, so a
///   scorer walks one contiguous run per genotype.
/// - **log multinomial coefficients** — `ln` of how many ways the genome's copies can
///   be ordered to spell the genotype: `ln 2` for a diploid heterozygote, `ln 1 = 0`
///   for a homozygote. The genotype prior multiplies each allele's frequency once per
///   copy and this is the count of orderings that product is missing
///   (`doc/devel/ng/spec/calling_priors.md` §3.1).
/// - **the homozygous lookup** — `Some(a)` where every copy is allele *a*, `None`
///   otherwise. **This is the one homozygous test in the caller**
///   (`doc/devel/ng/arch/calling_priors.md` §3.2): the inbreeding mixture's second
///   branch fires on it, and above diploidy what "homozygous" should mean is a known
///   open question deferred to a spec of its own
///   (`doc/devel/ng/spec/calling_priors.md` §3.3). Keeping it in one table is what
///   gives that spec one place to change — a second test written inline anywhere else
///   is where the two would silently diverge.
///
/// The fields are private and read through [`Self::view`] or the row lookups, because
/// the three tables are parallel — one entry, or one row, per genotype — and the
/// enumeration order is a contract with the VCF writer. A public field would let a
/// consumer hold two of them from different shapes.
#[derive(Clone, PartialEq, Debug)]
pub struct GenotypeTable {
    ploidy: Ploidy,
    allele_count: usize,
    genotype_count: usize,
    genotype_allele_counts: Vec<u32>,
    log_multinomial_coeffs: Vec<f64>,
    homozygous_allele_for: Vec<Option<AlleleId>>,
}

impl GenotypeTable {
    /// The table for one shape, from the cache when the shape is a common one and
    /// freshly built when it is not.
    ///
    /// Repeated calls for the same shape hand back the same `Arc` while the shape is
    /// within the cache bounds ([`MAX_CACHED_PLOIDY`], [`MAX_CACHED_ALLELE_COUNT`]);
    /// outside them every call builds again. Both answers hold the same values — see
    /// the constants for what the bound costs on either side of it.
    ///
    /// # Panics
    ///
    /// On `allele_count` of zero. A locus is called over at least its reference allele
    /// — `CandidateAlleles` cannot be built without one — so a table of no alleles is
    /// a shape derived from something other than a locus's allele table. It is
    /// refused here rather than returned empty, because an empty table has no
    /// homozygous-reference row and every consumer would read the absence as a locus
    /// with no candidate genotypes.
    ///
    /// On `allele_count` above [`MAX_ALLELE_COUNT`]. An [`AlleleId`] names exactly
    /// 65,536 alleles — ids 0 through 65,535 — so a 65,536-allele table is the widest
    /// one whose every row the homozygous lookup can still name; the 65,537th allele
    /// would need id 65,536, which an `AlleleId` does not have.
    ///
    /// On a shape whose table does not fit in a `usize`, whether because the genotype
    /// count itself does not or because that count times the allele count does not.
    /// Genotype count grows as `C(alleles + ploidy − 1, ploidy)`, so a wide locus at
    /// high ploidy is astronomically large long before it is merely big — 64 alleles
    /// at ploidy 40 is about 6.1e28, where a `usize` reaches about 1.8e19. Both
    /// halves are checked in one place, so what a caller sees names the shape rather
    /// than an arithmetic overflow inside an allocation.
    pub fn build(ploidy: Ploidy, allele_count: usize) -> Arc<Self> {
        let copies = usize::from(ploidy.get());
        if copies > MAX_CACHED_PLOIDY || allele_count > MAX_CACHED_ALLELE_COUNT {
            return Arc::new(Self::build_uncached(ploidy, allele_count));
        }
        let slot = copies * (MAX_CACHED_ALLELE_COUNT + 1) + allele_count;
        if let Some(table) = TABLE_CACHE.with(|cache| cache.borrow()[slot].clone()) {
            // The slot index and CACHE_SLOTS are written from the two bounds
            // separately, and an index that mixed them up would still land inside the
            // array — it would just hand one shape another's genotypes, which is a
            // wrong PL label rather than a crash.
            debug_assert_eq!(
                table.ploidy, ploidy,
                "cache slot {slot} holds a table of another ploidy"
            );
            debug_assert_eq!(
                table.allele_count, allele_count,
                "cache slot {slot} holds a table of another allele count"
            );
            return table;
        }
        // Built outside the borrow on purpose. Holding the cell across the build would
        // turn any future call that reached back into `build` — an instrumentation
        // hook, a richer shape type's `Display` — into a double-borrow panic, and the
        // release profile aborts on panic rather than unwinding one locus.
        let table = Arc::new(Self::build_uncached(ploidy, allele_count));
        TABLE_CACHE.with(|cache| {
            let mut slots = cache.borrow_mut();
            match &slots[slot] {
                Some(already) => Arc::clone(already),
                None => {
                    slots[slot] = Some(Arc::clone(&table));
                    table
                }
            }
        })
    }

    /// The build itself, without the cache — the same values [`Self::build`] returns.
    /// Separate so the cached and uncached paths cannot drift.
    fn build_uncached(ploidy: Ploidy, allele_count: usize) -> Self {
        assert!(
            allele_count > 0,
            "a locus is called over at least its reference allele, so a genotype table \
             of no alleles is a shape that came from somewhere other than a locus"
        );
        assert!(
            allele_count <= MAX_ALLELE_COUNT,
            "a genotype table names its alleles with AlleleId, which reaches 65,535: \
             {allele_count} alleles is past what it can name"
        );
        let genotype_count = count_genotypes(ploidy, allele_count)
            .and_then(|count| count.checked_mul(allele_count).map(|_| count))
            .unwrap_or_else(|| {
                panic!(
                    "the genotype table at ploidy {ploidy} over {allele_count} alleles \
                     does not fit in a usize, so it could never be built"
                )
            });
        let entries = genotype_count * allele_count;

        let genotype_allele_counts = enumerate_allele_counts(ploidy, allele_count, entries);
        // The count is a closed form and the rows come from a recursion: two
        // independent computations of one quantity, and a disagreement would be a
        // mis-sized likelihood table downstream rather than a crash here.
        assert_eq!(
            genotype_allele_counts.len(),
            entries,
            "the genotype count and the enumeration disagree at ploidy {ploidy} over \
             {allele_count} alleles"
        );

        let log_multinomial_coeffs: Vec<f64> = genotype_allele_counts
            .chunks_exact(allele_count)
            .map(|counts| log_multinomial_coefficient(ploidy, counts))
            .collect();

        let homozygous_allele_for: Vec<Option<AlleleId>> = genotype_allele_counts
            .chunks_exact(allele_count)
            .map(homozygous_allele)
            .collect();

        Self {
            ploidy,
            allele_count,
            genotype_count,
            genotype_allele_counts,
            log_multinomial_coeffs,
            homozygous_allele_for,
        }
    }

    /// The flat borrow the genotype prior and the read likelihood both take — the
    /// three tables and the shape they describe, in one argument
    /// (`doc/devel/ng/arch/calling_em_loop.md` §2; the prior's seam takes the three
    /// tables as separate flat slices, `doc/devel/ng/arch/calling_priors.md` §3.2, and
    /// reads them off this).
    ///
    /// `self` is destructured exhaustively rather than field by field, so that adding
    /// a seventh table to [`GenotypeTable`] — `nonzero_pairs` is the one the
    /// architecture names — fails to compile here instead of reaching consumers as a
    /// table the view silently does not carry.
    #[inline]
    pub fn view(&self) -> GenotypeTableView<'_> {
        let Self {
            ploidy,
            allele_count,
            genotype_count,
            genotype_allele_counts,
            log_multinomial_coeffs,
            homozygous_allele_for,
        } = self;
        GenotypeTableView {
            ploidy: *ploidy,
            allele_count: *allele_count,
            genotype_count: *genotype_count,
            genotype_allele_counts,
            log_multinomial_coeffs,
            homozygous_allele_for,
        }
    }

    /// How many copies of the genome a sample called against this table has.
    #[inline]
    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    /// How many alleles the locus is called over — the width of one allele-count row.
    #[inline]
    pub fn allele_count(&self) -> usize {
        self.allele_count
    }

    /// How many candidate genotypes the shape has, and so how many rows each of the
    /// three tables holds.
    #[inline]
    pub fn genotype_count(&self) -> usize {
        self.genotype_count
    }

    /// The whole allele-count table, `genotype_count × allele_count`, row-major.
    #[inline]
    pub fn genotype_allele_counts(&self) -> &[u32] {
        &self.genotype_allele_counts
    }

    /// One log multinomial coefficient per genotype, in genotype order.
    #[inline]
    pub fn log_multinomial_coeffs(&self) -> &[f64] {
        &self.log_multinomial_coeffs
    }

    /// The homozygous lookup, one entry per genotype — the caller's one homozygous
    /// test (see the type's documentation).
    #[inline]
    pub fn homozygous_alleles(&self) -> &[Option<AlleleId>] {
        &self.homozygous_allele_for
    }

    /// How many copies of each allele one genotype carries, or `None` if this table
    /// has no such genotype.
    ///
    /// Looked up rather than indexed, for the reason [`GenotypeIdx`] gives: an index
    /// minted against a different shape is a legal `u32` here, and indexing would
    /// either panic or hand back a real but wrong row.
    #[inline]
    pub fn allele_counts_of(&self, genotype: GenotypeIdx) -> Option<&[u32]> {
        self.view().allele_counts_of(genotype)
    }

    /// `ln` of the number of ways the genome's copies order to spell one genotype, or
    /// `None` if this table has no such genotype.
    #[inline]
    pub fn log_multinomial_coeff_of(&self, genotype: GenotypeIdx) -> Option<f64> {
        self.view().log_multinomial_coeff_of(genotype)
    }

    /// The allele one genotype is homozygous for — `Some(a)` when every copy is allele
    /// *a*, `None` when the genotype carries more than one allele, and `None` too when
    /// this table has no such genotype. **The two `None`s are not the same fact**, so
    /// a caller that must tell "not homozygous" from "no such row" reads
    /// [`Self::genotype_count`] first; the loop never has to, because its indices come
    /// from this table.
    #[inline]
    pub fn homozygous_allele_of(&self, genotype: GenotypeIdx) -> Option<AlleleId> {
        self.view().homozygous_allele_of(genotype)
    }
}

/// A borrow of one [`GenotypeTable`]'s three flat tables, handed across the calling
/// loop's two seams — the genotype prior fills a row of log priors from it, the read
/// likelihood a row of log likelihoods
/// (`doc/devel/ng/arch/calling_em_loop.md` §2, `doc/devel/ng/arch/calling_priors.md`
/// §3.2).
///
/// The fields are private and read through the accessors for the same reason the
/// table's are: the three slices are parallel and only [`GenotypeTable::view`] can
/// pair them. There is deliberately no public constructor, so a view always describes
/// a table that exists.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct GenotypeTableView<'a> {
    ploidy: Ploidy,
    allele_count: usize,
    genotype_count: usize,
    genotype_allele_counts: &'a [u32],
    log_multinomial_coeffs: &'a [f64],
    homozygous_allele_for: &'a [Option<AlleleId>],
}

impl<'a> GenotypeTableView<'a> {
    /// How many copies of the genome a sample called against this table has.
    #[inline]
    pub fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    /// How many alleles the locus is called over — the width of one allele-count row.
    #[inline]
    pub fn allele_count(&self) -> usize {
        self.allele_count
    }

    /// How many candidate genotypes the shape has.
    #[inline]
    pub fn genotype_count(&self) -> usize {
        self.genotype_count
    }

    /// The whole allele-count table, `genotype_count × allele_count`, row-major.
    #[inline]
    pub fn genotype_allele_counts(&self) -> &'a [u32] {
        self.genotype_allele_counts
    }

    /// One log multinomial coefficient per genotype, in genotype order.
    #[inline]
    pub fn log_multinomial_coeffs(&self) -> &'a [f64] {
        self.log_multinomial_coeffs
    }

    /// The homozygous lookup, one entry per genotype — the caller's one homozygous
    /// test (see [`GenotypeTable`]).
    #[inline]
    pub fn homozygous_alleles(&self) -> &'a [Option<AlleleId>] {
        self.homozygous_allele_for
    }

    /// The row one index names, or `None` if this shape has no such genotype. The one
    /// place an index becomes a row, so the three lookups below cannot come to
    /// disagree about which indices a table holds.
    #[inline]
    fn row_of(&self, genotype: GenotypeIdx) -> Option<usize> {
        let row = usize::try_from(genotype.get()).ok()?;
        (row < self.genotype_count).then_some(row)
    }

    /// How many copies of each allele one genotype carries, or `None` if this table
    /// has no such genotype.
    #[inline]
    pub fn allele_counts_of(&self, genotype: GenotypeIdx) -> Option<&'a [u32]> {
        let start = self.row_of(genotype)? * self.allele_count;
        self.genotype_allele_counts
            .get(start..start + self.allele_count)
    }

    /// `ln` of the number of ways the genome's copies order to spell one genotype, or
    /// `None` if this table has no such genotype.
    #[inline]
    pub fn log_multinomial_coeff_of(&self, genotype: GenotypeIdx) -> Option<f64> {
        self.log_multinomial_coeffs
            .get(self.row_of(genotype)?)
            .copied()
    }

    /// The allele one genotype is homozygous for; see
    /// [`GenotypeTable::homozygous_allele_of`] for what the two kinds of `None` mean.
    #[inline]
    pub fn homozygous_allele_of(&self, genotype: GenotypeIdx) -> Option<AlleleId> {
        self.homozygous_allele_for
            .get(self.row_of(genotype)?)
            .copied()
            .flatten()
    }
}

/// How many genotypes a shape has: `C(alleles + ploidy − 1, ploidy)`, the number of
/// multisets of `ploidy` alleles drawn from `allele_count` of them. `None` when that
/// does not fit in a `usize`.
///
/// Assumes a non-zero `allele_count`; [`GenotypeTable::build`] refuses zero before
/// this runs, and zero would come back from here as `None` — the right answer for the
/// wrong reason.
///
/// Computed by multiplying and dividing alternately, so every intermediate is itself a
/// binomial coefficient and the running value is as small as it can be — the naive
/// route through `(alleles + ploidy − 1)!` overflows at inputs this one handles
/// comfortably.
fn count_genotypes(ploidy: Ploidy, allele_count: usize) -> Option<usize> {
    let copies = u128::from(ploidy.get());
    let drawn_from = u128::try_from(allele_count)
        .ok()?
        .checked_add(copies)?
        .checked_sub(1)?;
    let mut count: u128 = 1;
    for step in 1..=copies {
        count = count.checked_mul(drawn_from.checked_sub(copies)?.checked_add(step)?)? / step;
    }
    usize::try_from(count).ok()
}

/// The allele-count table, `entries` numbers in all — one row of `allele_count` per
/// genotype — in the VCF genotype order, the order production's `genotype_order`
/// (`src/var_calling/per_group_merger.rs:522`) produces. Pinned by hand and as a law
/// here; pinned against that function in step C2.
fn enumerate_allele_counts(ploidy: Ploidy, allele_count: usize, entries: usize) -> Vec<u32> {
    let mut rows = Vec::with_capacity(entries);
    let mut counts = vec![0_u32; allele_count];
    push_genotypes_with_highest_allele_below(
        usize::from(ploidy.get()),
        allele_count,
        &mut counts,
        &mut rows,
    );
    rows
}

/// Emit every genotype of `copies_left` copies whose alleles are all below
/// `alleles_allowed`, appending one allele-count row per genotype.
///
/// **The order is what makes this the VCF order.** Each level picks the *highest*
/// allele the genotype carries and recurses on the copies below it, so genotypes come
/// out grouped by their highest allele, and within a group ordered by the same rule
/// one copy down. At ploidy 2 over three alleles that is `0/0`, `0/1`, `1/1`, `0/2`,
/// `1/2`, `2/2` — the order `PL` is written in. Building the order directly costs no
/// sort and, more to the point, leaves no comparator to get subtly wrong.
fn push_genotypes_with_highest_allele_below(
    copies_left: usize,
    alleles_allowed: usize,
    counts: &mut [u32],
    rows: &mut Vec<u32>,
) {
    if copies_left == 0 {
        rows.extend_from_slice(counts);
        return;
    }
    for highest in 0..alleles_allowed {
        counts[highest] += 1;
        push_genotypes_with_highest_allele_below(copies_left - 1, highest + 1, counts, rows);
        counts[highest] -= 1;
    }
}

/// `ln(n!)`, summed term by term.
///
/// The summation order is `ln 2 + ln 3 + … + ln n`, and it is kept because production
/// computes the same quantity the same way: the two are compared value for value, and
/// a different order would differ in the last bits. `n` never exceeds the ploidy, so
/// this is at most a handful of `ln` calls.
fn log_factorial(n: u32) -> f64 {
    let mut acc = 0.0_f64;
    for i in 2..=n {
        acc += f64::from(i).ln();
    }
    acc
}

/// `ln` of the multinomial coefficient `ploidy! / ∏ counts!` — how many ways the
/// genome's copies can be ordered to spell a genotype, in logs.
fn log_multinomial_coefficient(ploidy: Ploidy, counts: &[u32]) -> f64 {
    let mut log_coefficient = log_factorial(u32::from(ploidy.get()));
    for &copies in counts {
        log_coefficient -= log_factorial(copies);
    }
    log_coefficient
}

/// The allele every copy is, when there is one — `Some(a)` if allele *a* is the only
/// one with a non-zero count, `None` otherwise. See [`GenotypeTable`] for why this is
/// the caller's single homozygous test.
fn homozygous_allele(counts: &[u32]) -> Option<AlleleId> {
    let mut only_allele_seen: Option<AlleleId> = None;
    for (allele, &copies) in counts.iter().enumerate() {
        if copies == 0 {
            continue;
        }
        if only_allele_seen.is_some() {
            return None;
        }
        // PANIC-FREE: `allele` indexes a row of `counts`, whose width is the table's
        // `allele_count`, and the build refuses an `allele_count` above
        // `MAX_ALLELE_COUNT` — so the index is at most 65,535 and always fits a `u16`.
        only_allele_seen = Some(AlleleId(
            u16::try_from(allele)
                .expect("a genotype table refuses more alleles than AlleleId names"),
        ));
    }
    only_allele_seen
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("test ploidy is non-zero")
    }

    /// Rows of a table as tuples of allele ids, which is how a reader thinks of a
    /// genotype — `[2, 0]` is `0/0` and `[1, 1]` is `0/1`.
    fn genotypes_as_allele_lists(table: &GenotypeTable) -> Vec<Vec<u16>> {
        table
            .genotype_allele_counts()
            .chunks_exact(table.allele_count())
            .map(|counts| {
                counts
                    .iter()
                    .enumerate()
                    .flat_map(|(allele, &copies)| {
                        std::iter::repeat_n(allele as u16, copies as usize)
                    })
                    .collect()
            })
            .collect()
    }

    // ---- The enumeration ----

    /// The order `PL` and `GL` are written in. A table that enumerated differently
    /// would label every likelihood with the wrong genotype and nothing would crash.
    #[test]
    fn diploid_triallelic_genotypes_come_out_in_vcf_order() {
        let table = GenotypeTable::build(ploidy(2), 3);
        assert_eq!(
            genotypes_as_allele_lists(&table),
            vec![
                vec![0, 0],
                vec![0, 1],
                vec![1, 1],
                vec![0, 2],
                vec![1, 2],
                vec![2, 2],
            ]
        );
    }

    /// The same order as a **law that holds at every shape**, rather than at the two
    /// shapes written out by hand: consecutive rows, read as sorted allele lists, are
    /// strictly increasing when compared from the highest copy downwards. That is the
    /// comparator production's `genotype_order` sorts with
    /// (`src/var_calling/per_group_merger.rs`), and it is what pins the order at
    /// ploidy 3 and above over three or more alleles — where swapping two rows changes
    /// which one is homozygous-reference while every other test here still passes,
    /// because the three tables permute together.
    #[test]
    fn every_consecutive_pair_of_rows_is_in_vcf_order() {
        for copies in [1_u8, 2, 3, 4, 6, 8] {
            for allele_count in 1..=6_usize {
                let table = GenotypeTable::build(ploidy(copies), allele_count);
                for pair in genotypes_as_allele_lists(&table).windows(2) {
                    let (before, after) = (&pair[0], &pair[1]);
                    let order = before
                        .iter()
                        .rev()
                        .zip(after.iter().rev())
                        .find_map(|(here, there)| match here.cmp(there) {
                            std::cmp::Ordering::Equal => None,
                            other => Some(other),
                        })
                        .unwrap_or(std::cmp::Ordering::Equal);
                    assert_eq!(
                        order,
                        std::cmp::Ordering::Less,
                        "ploidy {copies}, {allele_count} alleles: {before:?} is not \
                         before {after:?} in the VCF order"
                    );
                }
            }
        }
    }

    #[test]
    fn diploid_biallelic_table_holds_the_three_genotypes_with_their_copy_counts() {
        let table = GenotypeTable::build(ploidy(2), 2);
        assert_eq!(table.genotype_count(), 3);
        assert_eq!(table.genotype_allele_counts(), &[2, 0, 1, 1, 0, 2]);
    }

    /// Ploidy above two is in this caller's stated range — a hexaploid or octoploid
    /// crop region — so the enumeration is pinned there too, and by hand.
    #[test]
    fn tetraploid_biallelic_table_runs_from_four_reference_copies_to_four_alternative() {
        let table = GenotypeTable::build(ploidy(4), 2);
        assert_eq!(
            genotypes_as_allele_lists(&table),
            vec![
                vec![0, 0, 0, 0],
                vec![0, 0, 0, 1],
                vec![0, 0, 1, 1],
                vec![0, 1, 1, 1],
                vec![1, 1, 1, 1],
            ]
        );
    }

    /// One sample with one copy of the genome: the genotypes are the alleles, and
    /// every one of them is homozygous by the definition the prior uses.
    #[test]
    fn a_haploid_locus_has_one_genotype_per_allele_and_all_of_them_are_homozygous() {
        let table = GenotypeTable::build(ploidy(1), 4);
        assert_eq!(table.genotype_count(), 4);
        assert_eq!(
            table.genotype_allele_counts(),
            &[1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
        );
        assert_eq!(
            table.homozygous_alleles(),
            &[
                Some(AlleleId(0)),
                Some(AlleleId(1)),
                Some(AlleleId(2)),
                Some(AlleleId(3)),
            ]
        );
    }

    /// Every genotype spends all of the sample's genome copies and no more: a row
    /// summing to anything but the ploidy is a genotype that is not one.
    #[test]
    fn every_genotype_row_accounts_for_exactly_the_ploidy_copies() {
        for copies in [1_u8, 2, 3, 4, 6, 8] {
            for allele_count in 1..=6_usize {
                let table = GenotypeTable::build(ploidy(copies), allele_count);
                for row in table.genotype_allele_counts().chunks_exact(allele_count) {
                    assert_eq!(
                        row.iter().sum::<u32>(),
                        u32::from(copies),
                        "ploidy {copies}, {allele_count} alleles, row {row:?}"
                    );
                }
            }
        }
    }

    /// The count is the number of multisets, `C(alleles + ploidy − 1, ploidy)`, and
    /// the enumeration has to produce exactly that many rows — no genotype emitted
    /// twice, none missed.
    #[test]
    fn the_table_holds_one_row_per_multiset_of_alleles() {
        // (ploidy, allele count, genotypes) — C(a + p - 1, p), by hand.
        let expected = [
            (1_u8, 1_usize, 1_usize),
            (1, 6, 6),
            (2, 1, 1),
            (2, 2, 3),
            (2, 3, 6),
            (2, 6, 21),
            (4, 2, 5),
            (4, 6, 126),
            (8, 16, 490_314),
        ];
        for (copies, allele_count, genotypes) in expected {
            let table = GenotypeTable::build(ploidy(copies), allele_count);
            assert_eq!(
                table.genotype_count(),
                genotypes,
                "ploidy {copies}, {allele_count} alleles"
            );
            assert_eq!(
                table.genotype_allele_counts().len(),
                genotypes * allele_count
            );
            assert_eq!(table.log_multinomial_coeffs().len(), genotypes);
            assert_eq!(table.homozygous_alleles().len(), genotypes);
        }
    }

    /// No genotype appears twice — the check the count alone cannot make, since a
    /// duplicate and a miss would cancel.
    #[test]
    fn no_genotype_is_enumerated_twice() {
        for copies in [1_u8, 2, 4, 8] {
            for allele_count in 1..=6_usize {
                let table = GenotypeTable::build(ploidy(copies), allele_count);
                let mut rows = genotypes_as_allele_lists(&table);
                let enumerated = rows.len();
                rows.sort();
                rows.dedup();
                assert_eq!(
                    rows.len(),
                    enumerated,
                    "ploidy {copies}, {allele_count} alleles: a genotype was emitted twice"
                );
            }
        }
    }

    // ---- The coefficients ----

    /// `ln 1`, `ln 2`, `ln 1` — one ordering for each homozygote, two for the
    /// heterozygote.
    #[test]
    fn diploid_biallelic_coefficients_are_one_two_one_in_logs() {
        let table = GenotypeTable::build(ploidy(2), 2);
        assert_eq!(
            table.log_multinomial_coeffs(),
            &[0.0, std::f64::consts::LN_2, 0.0]
        );
    }

    /// The coefficient counts orderings, so it must equal `ln` of the exact integer
    /// `ploidy! / ∏ counts!` — computed here the other way round, from integer
    /// factorials, so the check does not repeat the implementation's own summation.
    #[test]
    fn every_coefficient_is_the_log_of_the_exact_number_of_orderings() {
        fn factorial(n: u32) -> u128 {
            (1..=u128::from(n)).product::<u128>()
        }
        for copies in [1_u8, 2, 3, 4, 6, 8] {
            for allele_count in 1..=6_usize {
                let table = GenotypeTable::build(ploidy(copies), allele_count);
                for (row, counts) in table
                    .genotype_allele_counts()
                    .chunks_exact(allele_count)
                    .enumerate()
                {
                    let orderings = counts
                        .iter()
                        .fold(factorial(u32::from(copies)), |acc, &k| acc / factorial(k));
                    let expected = (orderings as f64).ln();
                    let got = table.log_multinomial_coeffs()[row];
                    assert!(
                        (got - expected).abs() < 1e-12,
                        "ploidy {copies}, {allele_count} alleles, row {row} {counts:?}: \
                         got {got}, exact ln({orderings}) is {expected}"
                    );
                }
            }
        }
    }

    // ---- The homozygous lookup ----

    /// Exactly the all-one-allele rows are homozygous, and each names its own allele.
    #[test]
    fn the_homozygous_lookup_names_an_allele_exactly_where_every_copy_is_that_allele() {
        for copies in [1_u8, 2, 3, 4, 8] {
            for allele_count in 1..=6_usize {
                let table = GenotypeTable::build(ploidy(copies), allele_count);
                for (row, counts) in table
                    .genotype_allele_counts()
                    .chunks_exact(allele_count)
                    .enumerate()
                {
                    let all_copies_at = counts
                        .iter()
                        .position(|&k| k == u32::from(copies))
                        .map(|allele| AlleleId(allele as u16));
                    assert_eq!(
                        table.homozygous_alleles()[row],
                        all_copies_at,
                        "ploidy {copies}, {allele_count} alleles, row {row} {counts:?}"
                    );
                }
                // One homozygote per allele, whatever the ploidy.
                assert_eq!(
                    table
                        .homozygous_alleles()
                        .iter()
                        .filter(|entry| entry.is_some())
                        .count(),
                    allele_count
                );
            }
        }
    }

    /// The diploid triallelic case spelled out, because it is the one a reader can
    /// check against the VCF order by eye: `0/0`, `0/1`, `1/1`, `0/2`, `1/2`, `2/2`.
    #[test]
    fn diploid_triallelic_homozygous_lookup_marks_the_first_third_and_last_genotypes() {
        let table = GenotypeTable::build(ploidy(2), 3);
        assert_eq!(
            table.homozygous_alleles(),
            &[
                Some(AlleleId(0)),
                None,
                Some(AlleleId(1)),
                None,
                None,
                Some(AlleleId(2)),
            ]
        );
    }

    // ---- The row lookups ----

    #[test]
    fn a_row_lookup_returns_that_genotypes_counts_coefficient_and_homozygous_allele() {
        let table = GenotypeTable::build(ploidy(2), 3);
        // Row 2 of the VCF order is 1/1.
        assert_eq!(table.allele_counts_of(GenotypeIdx(2)), Some(&[0, 2, 0][..]));
        assert_eq!(table.log_multinomial_coeff_of(GenotypeIdx(2)), Some(0.0));
        assert_eq!(
            table.homozygous_allele_of(GenotypeIdx(2)),
            Some(AlleleId(1))
        );
        // Row 4 is 1/2 — carries two alleles, so no homozygous allele.
        assert_eq!(table.allele_counts_of(GenotypeIdx(4)), Some(&[0, 1, 1][..]));
        assert_eq!(
            table.log_multinomial_coeff_of(GenotypeIdx(4)),
            Some(std::f64::consts::LN_2)
        );
        assert_eq!(table.homozygous_allele_of(GenotypeIdx(4)), None);
    }

    /// An index minted against a wider shape is a legal `u32` here, and the lookup
    /// must say so rather than hand back a real but wrong row.
    #[test]
    fn a_row_lookup_past_the_end_of_the_table_finds_nothing() {
        let table = GenotypeTable::build(ploidy(2), 2);
        assert_eq!(table.genotype_count(), 3);
        assert_eq!(table.allele_counts_of(GenotypeIdx(3)), None);
        assert_eq!(table.log_multinomial_coeff_of(GenotypeIdx(3)), None);
        assert_eq!(table.homozygous_allele_of(GenotypeIdx(3)), None);
        assert_eq!(table.allele_counts_of(GenotypeIdx(u32::MAX)), None);
        assert_eq!(table.log_multinomial_coeff_of(GenotypeIdx(u32::MAX)), None);
        assert_eq!(table.homozygous_allele_of(GenotypeIdx(u32::MAX)), None);
    }

    // ---- The view ----

    #[test]
    fn the_view_carries_the_same_shape_and_the_same_three_tables() {
        let table = GenotypeTable::build(ploidy(4), 3);
        let view = table.view();
        assert_eq!(view.ploidy(), table.ploidy());
        assert_eq!(view.allele_count(), table.allele_count());
        assert_eq!(view.genotype_count(), table.genotype_count());
        assert_eq!(
            view.genotype_allele_counts(),
            table.genotype_allele_counts()
        );
        assert_eq!(
            view.log_multinomial_coeffs(),
            table.log_multinomial_coeffs()
        );
        assert_eq!(view.homozygous_alleles(), table.homozygous_alleles());
        // The three tables are parallel: one row, one entry, one entry per genotype.
        assert_eq!(
            view.genotype_allele_counts().len(),
            view.genotype_count() * view.allele_count()
        );
        assert_eq!(view.log_multinomial_coeffs().len(), view.genotype_count());
        assert_eq!(view.homozygous_alleles().len(), view.genotype_count());
    }

    // ---- The cache ----

    /// The point of the cache: a run meets one shape at millions of loci and builds it
    /// once per thread.
    #[test]
    fn two_builds_of_one_cached_shape_return_the_same_table() {
        let first = GenotypeTable::build(ploidy(2), 2);
        let second = GenotypeTable::build(ploidy(2), 2);
        assert!(Arc::ptr_eq(&first, &second));
    }

    #[test]
    fn different_shapes_get_different_tables() {
        let diploid_biallelic = GenotypeTable::build(ploidy(2), 2);
        let diploid_triallelic = GenotypeTable::build(ploidy(2), 3);
        let tetraploid_biallelic = GenotypeTable::build(ploidy(4), 2);
        assert!(!Arc::ptr_eq(&diploid_biallelic, &diploid_triallelic));
        assert!(!Arc::ptr_eq(&diploid_biallelic, &tetraploid_biallelic));
        assert_ne!(
            diploid_biallelic.genotype_count(),
            diploid_triallelic.genotype_count()
        );
    }

    /// **Every shape the cache keeps gets a table of its own shape back.** The slot
    /// index and the array's size are written from the two bounds separately, and an
    /// index that mixed them up would still land inside the array — ploidy 1 over 16
    /// alleles and ploidy 2 over 7 would share a slot — so one locus would be scored
    /// against another shape's genotypes, with nothing crashing. Walking the whole
    /// grid in one thread is what makes such a collision visible.
    #[test]
    fn every_cached_shape_gets_back_a_table_of_its_own_shape() {
        for copies in 1..=MAX_CACHED_PLOIDY {
            for allele_count in 1..=MAX_CACHED_ALLELE_COUNT {
                let table = GenotypeTable::build(ploidy(copies as u8), allele_count);
                assert_eq!(
                    usize::from(table.ploidy().get()),
                    copies,
                    "asked for ploidy {copies} over {allele_count} alleles"
                );
                assert_eq!(
                    table.allele_count(),
                    allele_count,
                    "asked for ploidy {copies} over {allele_count} alleles"
                );
            }
        }
    }

    /// Past the cache bounds the values are still right; only the sharing is lost.
    /// Pinned in both directions, because the bound is two comparisons and either
    /// could be dropped without any other test noticing.
    #[test]
    fn a_shape_past_the_cache_bounds_is_rebuilt_each_time_with_the_same_values() {
        for (copies, allele_count) in [(9_u8, 2_usize), (2, 17)] {
            let first = GenotypeTable::build(ploidy(copies), allele_count);
            let second = GenotypeTable::build(ploidy(copies), allele_count);
            assert!(
                !Arc::ptr_eq(&first, &second),
                "ploidy {copies} over {allele_count} alleles is outside the cache bounds"
            );
            assert_eq!(first, second);
        }
    }

    /// The uncached branch is the one every polyploid or wide locus takes, and
    /// comparing two of its runs with each other cannot catch a wrong value — the
    /// builder is a pure function, so it agrees with itself whatever it computes.
    /// These are the values themselves.
    ///
    /// **Row 1 is here because row 0 cannot fail.** At ploidy 9 over two alleles row 0
    /// is `[9, 0]`, whose coefficient is `ln 9! − ln 9!` — zero under the right formula
    /// and zero under a `log_factorial` capped at 8, or at 4, or at 2. Row 1 is
    /// `[8, 1]`, worth `ln 9`, and it is the row that separates them.
    #[test]
    fn a_shape_past_the_cache_bounds_holds_the_right_values() {
        let deep = GenotypeTable::build(ploidy(9), 2);
        assert_eq!(deep.genotype_count(), 10);
        assert_eq!(deep.genotype_allele_counts()[..4], [9, 0, 8, 1]);
        assert_eq!(deep.log_multinomial_coeffs()[0], 0.0);
        // `ln 9`, to within a few units in the last place: the coefficient is summed as
        // `ln 9! − ln 8!` rather than computed as `ln 9`, so the two agree in value and
        // not in every bit. The exact bits are pinned against production in
        // `genotype_table_parity`; what this row is here to catch is a coefficient off
        // by whole nats.
        assert!(
            (deep.log_multinomial_coeffs()[1] - 9.0_f64.ln()).abs() < 1e-12,
            "ploidy 9, row 1 [8, 1] should be ln 9, got {}",
            deep.log_multinomial_coeffs()[1]
        );
        assert_eq!(deep.homozygous_alleles()[0], Some(AlleleId(0)));
        assert_eq!(deep.homozygous_alleles()[1], None);
        assert_eq!(deep.homozygous_alleles()[9], Some(AlleleId(1)));

        let wide = GenotypeTable::build(ploidy(2), 17);
        assert_eq!(wide.genotype_count(), 153);
        assert_eq!(wide.log_multinomial_coeffs()[1], std::f64::consts::LN_2);
        assert_eq!(wide.homozygous_alleles()[152], Some(AlleleId(16)));
    }

    /// The cache is a store, not a second implementation: what it hands back has to be
    /// what the builder produces for the same shape.
    #[test]
    fn build_returns_the_same_values_from_the_cache_and_from_a_fresh_build() {
        for (copies, allele_count) in [(1_u8, 1_usize), (2, 2), (4, 6), (2, 16)] {
            let cached = GenotypeTable::build(ploidy(copies), allele_count);
            let fresh = GenotypeTable::build_uncached(ploidy(copies), allele_count);
            assert_eq!(*cached, fresh, "ploidy {copies}, {allele_count} alleles");
        }
    }

    /// The last shape the cache does keep, on both axes — the boundary the
    /// past-the-bounds test is one step beyond.
    #[test]
    fn the_widest_and_deepest_cached_shapes_are_still_shared() {
        for (copies, allele_count) in [(8_u8, 2_usize), (2, 16)] {
            let first = GenotypeTable::build(ploidy(copies), allele_count);
            let second = GenotypeTable::build(ploidy(copies), allele_count);
            assert!(
                Arc::ptr_eq(&first, &second),
                "ploidy {copies} over {allele_count} alleles is inside the cache bounds"
            );
        }
    }

    /// Each thread builds its own, so no lock sits on the loop's path — and a table
    /// cached by one thread is not handed to another. The negative pointer comparison
    /// is sound because `here` is still alive while `there` is built, so the two
    /// allocations cannot coincide.
    #[test]
    fn each_thread_keeps_its_own_cached_tables() {
        let here = GenotypeTable::build(ploidy(2), 4);
        let there = std::thread::spawn(|| GenotypeTable::build(ploidy(2), 4))
            .join()
            .expect("the builder thread does not panic");
        assert!(!Arc::ptr_eq(&here, &there));
        assert_eq!(here, there);
    }

    // ---- Refusals ----

    #[test]
    #[should_panic(expected = "called over at least its reference allele")]
    fn a_table_of_no_alleles_is_refused() {
        let _ = GenotypeTable::build(ploidy(2), 0);
    }

    #[test]
    #[should_panic(expected = "past what it can name")]
    fn a_table_wider_than_an_allele_id_can_name_is_refused() {
        let _ = GenotypeTable::build(ploidy(1), MAX_ALLELE_COUNT + 1);
    }

    #[test]
    #[should_panic(expected = "does not fit in a usize")]
    fn a_shape_with_more_genotypes_than_a_usize_can_count_is_refused() {
        // C(64 + 40 - 1, 40) is about 6.1e28; a usize reaches 2^64, about 1.8e19.
        // This shape exits through the final `usize::try_from`.
        let _ = GenotypeTable::build(ploidy(40), 64);
    }

    /// The other way the count can fail: the running binomial itself overflows before
    /// the final narrowing is ever reached.
    #[test]
    #[should_panic(expected = "does not fit in a usize")]
    fn a_shape_that_overflows_the_running_binomial_is_refused() {
        let _ = GenotypeTable::build(ploidy(255), MAX_ALLELE_COUNT);
    }

    /// A shape whose genotype *count* fits a `usize` while the flat table does not.
    /// At ploidy 4 over 65,536 alleles the count is 768,684,707,117,285,376 — which
    /// fits — and that many rows of 65,536 entries does not.
    #[test]
    #[should_panic(expected = "does not fit in a usize")]
    fn a_shape_whose_flat_table_outgrows_a_usize_is_refused() {
        let _ = GenotypeTable::build(ploidy(4), MAX_ALLELE_COUNT);
    }
}
