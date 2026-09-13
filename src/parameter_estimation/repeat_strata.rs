//! **How a repeat tract's stratum is named** — a motif period, the reference's repeat count, and
//! the library and ploidy that make one fitted set of slippage numbers.
//!
//! A *stratum* is one group of repeat tracts that gets its own fitted slippage parameters: every
//! tract sharing a motif period and a reference repeat count. Slippage depends on repeat count
//! more than on anything else, which is what makes that the grouping.
//!
//! **These three types are vocabulary rather than machinery.** They name what a fit is about; they
//! do no fitting. The fitting lives in [`joint::ssr_fit`](super::joint::ssr_fit), which reads a
//! cohort's censuses, and the calling step reads the results back under these keys
//! (`calling::inference::repeat_tract_parameters`).
//!
//! **There is a second `Stratum` in this tree and the difference is deliberate.**
//! [`joint::census::Stratum`](super::joint::census::Stratum) names a stratum with plain numbers —
//! a period in bases and a repeat count as a `u64` — because that is what a census record can
//! hold, and `run::census_fit` converts at the seam, dropping any stratum the checked types here
//! reject rather than coercing it. Two spellings of one concept is a cost; a census record that
//! carried checked types would be a format that could not be read back by a build whose ranges had
//! moved, which is worse.
//!
//! *Lifted out of `parameter_estimation::ssr` when the whole-genome histogram route was removed
//! (`impl_plan/remove_histogram_route.md`), which is where they used to live beside the
//! per-sample STR fit that is gone.*

use std::fmt;

use crate::types::{Ploidy, ReadGroupId, SsrPeriod};

/// How many whole motif copies a tract holds.
///
/// **The reference tract's count, never the sample's.** It is a pure function of the reference,
/// which is what makes every sample file a locus under the same stratum and so lets a cohort
/// compare one sample's slippage with another's (`arch/parameter_prepass_ssr.md` §2.1). A sample
/// whose alleles differ from the reference does not move between strata; its reads land at an
/// offset from the reference length instead.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct RepeatCount(pub u32);

impl RepeatCount {
    #[inline]
    #[must_use]
    pub fn get(self) -> u32 {
        self.0
    }
}

impl fmt::Display for RepeatCount {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// One group of loci that gets its own fitted slippage parameters: a motif period and a reference
/// repeat count (`spec/parameter_prepass_ssr.md` §4).
///
/// **Ordered by `(period, repeats)`**, so walking a map of strata visits each period's repeat
/// counts in ascending order — which is what the monotonicity rule needs: slippage genuinely rises
/// with repeat count, so a fitted sequence that dips in the middle is reporting noise in one
/// stratum rather than a fact about repeats, and finding that means comparing each stratum with
/// the one before it (§4.3).
///
/// **The reference's own catalog already speaks this concept and spells it differently**:
/// [`repeat_catalog::strata`](crate::repeat_catalog::strata) keys its counts and its
/// per-stratum sample on a raw `(u8, u64)` pair, so a driver that asks the catalog how many loci a
/// stratum holds converts at that seam. The two orderings agree — both are period first, then
/// repeat count.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Stratum {
    pub period: SsrPeriod,
    pub repeats: RepeatCount,
}

impl Stratum {
    #[must_use]
    pub fn new(period: SsrPeriod, repeats: RepeatCount) -> Self {
        Self { period, repeats }
    }
}

impl fmt::Display for Stratum {
    /// "period 2, 6 repeats" — the words the emitted summary and every error message use, so
    /// neither has to spell the pair out again.
    ///
    /// Destructured rather than read field by field, so that a third field added to [`Stratum`] is
    /// a compile error here rather than a rendering that silently names two thirds of the key.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self { period, repeats } = self;
        write!(f, "period {period}, {repeats} repeats")
    }
}

/// What one fitted set of slippage numbers is *about*: a library, a stratum, and how many genome
/// copies its loci sit on.
///
/// **The ploidy is in the key and not merely carried alongside it**, and that is the one part of
/// this type worth arguing. Evidence at a tract is the same object whatever the ploidy — a count
/// of reads at each offset, which knows nothing about how many chromosomes produced them — so
/// pooling a haploid locus with a diploid one is invisible while the loci are being counted. It
/// becomes wrong at the fit: the fit scores each entry against the genotypes of *one* ploidy, and
/// a pooled table has no ploidy that is true of all of it.
///
/// **It costs nothing on today's runs.** Every genome region is currently declared to have the
/// same ploidy, so every key carries the same value and the tables are exactly those a two-part
/// key would have built. What it buys is that the first sex chromosome or mixed-ploidy genome to
/// arrive splits the tables instead of silently merging them.
///
/// The field order is the order the fits walk it in: read group, then stratum, then ploidy.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct StratumKey {
    /// Which library. Slippage is a property of the chemistry, so each library is fitted
    /// separately.
    pub read_group: ReadGroupId,
    /// The motif period and the reference repeat count.
    pub stratum: Stratum,
    /// How many genome copies these loci sit on — the set of genotypes the fit scores each of
    /// this stratum's entries against.
    pub ploidy: Ploidy,
}

impl fmt::Display for StratumKey {
    /// For the error messages and the summary, which name a fit by all three of these.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "read group {}, {}, ploidy {}",
            self.read_group, self.stratum, self.ploidy
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn stratum(period: usize, repeats: u32) -> Stratum {
        Stratum::new(
            SsrPeriod::try_new(period).expect("a positive period"),
            RepeatCount(repeats),
        )
    }

    /// The ordering is what the monotonicity rule walks, so it is asserted rather than assumed:
    /// period first, then repeat count.
    #[test]
    fn strata_order_by_period_then_repeat_count() {
        let mut all = vec![stratum(2, 6), stratum(1, 30), stratum(2, 5), stratum(1, 8)];
        all.sort();
        assert_eq!(
            all,
            vec![stratum(1, 8), stratum(1, 30), stratum(2, 5), stratum(2, 6)]
        );
    }

    /// The words every error message and the run summary use.
    #[test]
    fn a_stratum_renders_as_the_words_the_summary_uses() {
        assert_eq!(stratum(2, 6).to_string(), "period 2, 6 repeats");
        assert_eq!(
            StratumKey {
                read_group: ReadGroupId(3),
                stratum: stratum(1, 12),
                ploidy: Ploidy::try_new(2).expect("two is a ploidy"),
            }
            .to_string(),
            "read group 3, period 1, 12 repeats, ploidy 2"
        );
    }
}
