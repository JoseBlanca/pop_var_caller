//! The STR path: how often a read at a repeat tract shows a length its allele does not
//! have, and the four numbers that describe it.
//!
//! A read at an ordinary site can only be misread — one base for another. A read at a
//! repeat tract can also **slip**, showing a whole motif copy more or fewer than the
//! allele it was drawn from, and no per-base error rate has anywhere to put that. So this
//! path carries the generic path's noise model and adds slippage to it: how often a read
//! slips at all, which way it slips, how far it slips when it does, and the per-base
//! substitution rate that changes a tract's composition at fixed length.
//!
//! **A stratum is a motif period and a reference repeat count** — one group of loci that
//! gets its own fitted numbers, because how much a tract slips depends on how many copies
//! it holds: 9 reads in 10,000 below four repeats against 2 in 100 at six or more, a
//! twenty-two-fold spread inside one dataset (`spec/parameter_prepass_ssr.md` §4). **The
//! fits are filed per `(read group, stratum)`**, the read group being the other half of the
//! key rather than part of the stratum, because slippage is a property of the chemistry
//! rather than of the individual.
//!
//! **Nothing here calls a genotype first.** The genotype is summed over rather than
//! chosen, which is what separates this from production's stutter pre-pass: that one pools
//! reads from loci that passed a confident-genotype gate, and measured against HG002's
//! assembly truth it reports gains as *more* common than losses where the truth is 3.4
//! times the other way (`spec/parameter_prepass.md` §2.2).
//!
//! Design: `doc/devel/ng/spec/parameter_prepass_ssr.md` (what is fitted and why),
//! `doc/devel/ng/arch/parameter_prepass_ssr.md` (types and interfaces), on the shared
//! framing of `doc/devel/ng/spec/parameter_prepass.md`.
//!
//! Five files, split so that the shaping of data and the mathematics on it never live
//! together:
//!
//! - this one — the vocabulary a stratum is keyed on, the widths the fit works inside, and
//!   (from Milestone A5) what a fit emits;
//! - [`offset_bucket`] — where a read's tract length is recorded, and the only place that
//!   can build one of those buckets;
//! - [`locus_offsets`] — one STR locus reduced to one table entry;
//! - [`stratum_table`] — the sparse table of locus shapes, per stratum;
//! - [`slippage`] — the noise model and the search that fits it.
//!
//! **No trait over the accumulator.** Nothing generic drives it and the walk knows which
//! object it is filling; `fitting/` is the one genuine swappable seam in step 4.

pub mod locus_offsets;
pub mod offset_bucket;
pub mod slippage;
pub mod stratum_table;

use std::fmt;

use crate::ng::types::SsrPeriod;

pub use offset_bucket::{OFFSET_BUCKETS, OFFSET_HALF_RANGE, OffsetBucket, bucket_of};

/// How many whole motif copies a tract holds.
///
/// **The reference tract's count, never the sample's.** It is a pure function of the
/// reference, which is what makes every sample file a locus under the same stratum and so
/// lets a cohort compare one sample's stutter with another's
/// (`arch/parameter_prepass_ssr.md` §2.1). A sample whose alleles differ from the
/// reference does not move between strata; its reads land at an offset instead
/// ([`WholeRepeatOffset`]).
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

/// One group of loci that gets its own fitted stutter parameters: a motif period and a
/// reference repeat count (`spec/parameter_prepass_ssr.md` §4).
///
/// **Ordered by `(period, repeats)`**, so walking a map of strata visits each period's
/// repeat counts in ascending order — which is what the monotonicity rule needs: slippage
/// genuinely rises with repeat count, so a fitted sequence that dips in the middle is
/// reporting noise in one stratum rather than a fact about repeats, and finding that means
/// comparing each stratum with the one before it (§4.3).
///
/// **The reference's own catalog already speaks this concept and spells it differently**:
/// [`repeat_catalog::strata`](crate::ng::repeat_catalog::strata) keys its counts and its
/// per-stratum sample on a raw `(u8, u64)` pair, so a driver that asks the catalog how many
/// loci a stratum holds converts at that seam. The two orderings agree — both are period
/// first, then repeat count — and nothing yet observes that they do, because nothing
/// converts between them until Milestone C.
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
    /// "period 2, 6 repeats" — the words the emitted summary and every error message use,
    /// so neither has to spell the pair out again.
    ///
    /// Destructured rather than read field by field, so that a third field added to
    /// `Stratum` is a compile error here rather than a rendering that silently names two
    /// thirds of the key.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let Self { period, repeats } = self;
        write!(f, "period {period}, {repeats} repeats")
    }
}

/// How far a read's tract sits from the **reference** tract length, in whole motif copies.
/// Negative is shorter than the reference, positive is longer.
///
/// **The origin is the reference tract's length and not each locus's most common observed
/// length**, which is what an earlier draft of the design used. Measured, the modal origin
/// returns a slippage level 50% to 408% high and destroys the direction asymmetry the model
/// exists to carry — a split of 0.48 where the truth is 0.17, a 1.1-fold imbalance where
/// the truth is 4.9-fold (`spec/parameter_prepass_ssr.md` §4.1). The reason is that
/// centring on the mode makes the origin a function of the reads, so a fit treating it as a
/// property of the locus is answering a question about a quantity that moved.
///
/// Unconstrained: any offset is a legal thing to observe. What is bounded is the bucket it
/// is *recorded* in ([`OffsetBucket`]) and the offsets a fitted allele may sit at
/// ([`ALLELE_OFFSET_LIMIT`]).
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct WholeRepeatOffset(pub i8);

impl WholeRepeatOffset {
    #[inline]
    #[must_use]
    pub fn get(self) -> i8 {
        self.0
    }
}

impl fmt::Display for WholeRepeatOffset {
    /// Signed, always — `+2` and `-2` are different observations and a bare `2` in a
    /// report would not say which.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:+}", self.0)
    }
}

/// How far from the reference tract length the fit may place an allele — the support of a
/// stratum's genotype frequencies, and **a different width from the recorded range, wider
/// than it**.
///
/// It is what lets a saturating end bucket be explained by a distant *allele* rather than
/// by a far *slip*, which is why this width and not [`OFFSET_HALF_RANGE`] decides the
/// answer. Too narrow and the fit charges a real long allele to slippage; too wide and a
/// thin stratum is fitting frequencies for lengths no locus carries, at `A(A+1)/2` numbers
/// for `A` lengths.
///
/// **Six comes from the measured distribution of allele lengths and is a threshold to
/// clear rather than a number to tune.** A locus's most common observed length against the
/// reference, on HG002 at 300×: 88.9% sit exactly at the reference, ±4 holds 99%, ±12 holds
/// 99.9%, ±19 holds 99.99%; tomato is tighter, 95.7% at zero with ±1 holding 99%
/// (`spec/parameter_prepass_ssr.md` §8, item 1). What a locus outside the support costs is
/// nothing and then everything: leaving 2.5% of loci out costs 0.1% of the slippage level,
/// 7.9% costs 2.5%, and **19.3% costs +499% with the direction asymmetry destroyed**. Six
/// leaves about one human locus in 200 outside — a fifth of the way to the row that is
/// already free.
///
/// **It clips at the low end, so the number of allele lengths is a per-stratum quantity**:
/// an allele cannot be shorter than nothing, so a stratum at 4 repeats reaches only −4.
/// [`allele_support`] is where that happens.
pub const ALLELE_OFFSET_LIMIT: i8 = 6;

/// The relation the design turns on, held by the compiler rather than by prose: the fit may
/// place an allele **further** from the reference than an entry records offsets. That is
/// what lets a saturating end bucket be attributed to a distant allele instead of to a far
/// slip; invert it and every locus whose allele sits past the recorded range has its reads
/// explained the only way left, as slippage.
const _: () = assert!(
    ALLELE_OFFSET_LIMIT > OFFSET_HALF_RANGE,
    "the fit must be able to place an allele beyond the offsets an entry records"
);

/// Reads entered from one locus. A deeper locus is entered from a **random subsample** of
/// its reads down to this, seeded from the locus's position so a region-sharded walk and a
/// single-threaded one keep the same reads and merging stays exact.
///
/// **Not the memory knob it looks like.** Measured on HG002 at 300× over the GIAB
/// tandem-repeat set, the uncapped table is 0.43 entries a locus — deep data deduplicates,
/// because most loci at a clean tract are "every read at the reference length" and what
/// separates two entries is mostly their depth (`spec/parameter_prepass_ssr.md` §4.1).
///
/// **Nor is it a correctness limit**: the scoring rule is exactly unbiased at every depth
/// tried, to 45 reads a locus. What it does decide is the width of an entry's counters — a
/// `u8` bucket count wraps silently above 255 — so the cap and that width are one decision.
/// Twelve is a low starting value whose only cost is the precision of the reads it drops.
///
/// Consumed from Milestone C, where a locus first becomes an entry.
pub const MAX_LOCUS_READS: u32 = 12;

/// Above this share of the reads that differ from the reference tract length, a stratum is
/// one this noise model does not describe: what moved the reads is ordinary indel rather
/// than repeat slippage, and a slippage rate fitted there is mostly mis-modelled indel
/// however many loci stood behind it (`spec/parameter_prepass_ssr.md` §5).
///
/// **One in ten, and the bands it separates are 0.9% against 33.8% and 58.5%** — a factor
/// of three either way, so a stratum crossing it is unambiguous. *Soft*: three bands of one
/// dataset, and the number moves if the per-stratum distribution turns out to be continuous
/// rather than the two clumps that table suggests.
///
/// Consumed from Milestone B, where the table starts counting the two denominators.
pub const GUARD_SHARE_LIMIT: f64 = 0.10;

/// The allele lengths a stratum's fit may place mass on, as offsets from the reference
/// tract length, in ascending order.
///
/// `±ALLELE_OFFSET_LIMIT` around the reference, **clipped at the low end because an allele
/// cannot be shorter than nothing**: a tract of 4 reference copies reaches −4, so it has 11
/// lengths where a tract of 6 or more has the full 13. The genotype count follows —
/// `A(A+1)/2` unordered pairs, so 66 against 91 — which is why the fit is handed a stratum's
/// support rather than assuming one (`arch/parameter_prepass_ssr.md` §2.1).
#[must_use]
pub fn allele_support(reference_repeats: RepeatCount) -> Vec<WholeRepeatOffset> {
    // How far below the reference an allele could reach at this stratum, before the fit's
    // own limit applies: a tract of `n` copies can lose at most `n` of them.
    let reach_down = i8::try_from(reference_repeats.get()).unwrap_or(i8::MAX);
    // Parenthesised deliberately: `-ALLELE_OFFSET_LIMIT.min(reach_down)` parses the same way
    // but reads as `(-6).min(reach_down)`, which at 3 reference repeats gives −6 — a support
    // reaching three copies below an empty tract, silently.
    let lowest = -(reach_down.min(ALLELE_OFFSET_LIMIT));
    (lowest..=ALLELE_OFFSET_LIMIT)
        .map(WholeRepeatOffset)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn period(bases: u8) -> SsrPeriod {
        SsrPeriod::try_new(usize::from(bases)).expect("a period inside the STR scope")
    }

    fn stratum(bases: u8, repeats: u32) -> Stratum {
        Stratum::new(period(bases), RepeatCount(repeats))
    }

    /// The monotonicity rule walks a period's strata in repeat-count order and compares each
    /// fit with the one before it, so this ordering is what puts neighbours next to each
    /// other. Sorting by repeat count first would interleave the periods and compare a
    /// dinucleotide's fit against a hexamer's.
    #[test]
    fn strata_sort_by_period_and_then_by_repeat_count() {
        let mut strata = vec![
            stratum(2, 12),
            stratum(1, 20),
            stratum(2, 6),
            stratum(1, 8),
            stratum(6, 4),
        ];
        strata.sort();

        assert_eq!(
            strata,
            vec![
                stratum(1, 8),
                stratum(1, 20),
                stratum(2, 6),
                stratum(2, 12),
                stratum(6, 4),
            ]
        );
    }

    /// A stratum names itself in the words the summary and the error messages use, so
    /// neither has to spell the pair out again — and a reader of a log sees the pair rather
    /// than a tuple.
    #[test]
    fn a_stratum_names_its_period_and_its_repeat_count() {
        assert_eq!(stratum(2, 6).to_string(), "period 2, 6 repeats");
    }

    /// **A zero period is unrepresentable here rather than rejected here**: `Stratum::new`
    /// takes an already-checked [`SsrPeriod`], and the rejection itself is
    /// `ng::types::tests::ssr_period_accepts_exactly_the_str_scope`. What this module owns is
    /// that a stratum carries back the pair it was built from — a swap between the two would
    /// file every locus under a stratum no model was fitted at.
    #[test]
    fn a_stratum_carries_back_the_period_and_repeat_count_it_was_built_from() {
        let built = stratum(2, 6);
        assert_eq!(built.period.get(), 2);
        assert_eq!(built.repeats.get(), 6);
    }

    /// An offset is read for its **direction** first — a report column of them is how one
    /// tells a stratum losing repeats from one gaining them — so the sign is explicit on
    /// gains and at the origin, not only on losses.
    #[test]
    fn an_offset_renders_with_its_sign_in_both_directions() {
        assert_eq!(WholeRepeatOffset(2).to_string(), "+2");
        assert_eq!(WholeRepeatOffset(-2).to_string(), "-2");
        assert_eq!(WholeRepeatOffset(0).to_string(), "+0");
    }

    /// **The support is every offset from the clip to the limit, with nothing missing and
    /// nothing repeated.** Asserting its length and its two ends is not enough: a support
    /// holding a duplicate and a hole has the right cardinality, so the `A(A+1)/2` genotype
    /// count still comes out right while the fit never considers one allele length and
    /// double-counts another.
    ///
    /// The five strata are the clip's whole range: none at 0 copies, biting at 3 and at 5,
    /// exactly released at 6, and far above it at 20.
    #[test]
    fn the_allele_support_is_every_offset_from_the_clip_to_the_limit() {
        for (repeats, lowest) in [(0u32, 0i8), (3, -3), (5, -5), (6, -6), (20, -6)] {
            assert_eq!(
                allele_support(RepeatCount(repeats)),
                (lowest..=ALLELE_OFFSET_LIMIT)
                    .map(WholeRepeatOffset)
                    .collect::<Vec<_>>(),
                "at {repeats} reference repeats"
            );
        }
    }

    /// The lengths behind the genotype count the design reasons about: 10, 13 and 13 allele
    /// lengths give 55, 91 and 91 unordered pairs.
    #[test]
    fn the_allele_support_clips_below_the_reference_but_not_above() {
        let short = allele_support(RepeatCount(3));
        assert_eq!(short.len(), 10);
        assert_eq!(short.first(), Some(&WholeRepeatOffset(-3)));
        assert_eq!(short.last(), Some(&WholeRepeatOffset(ALLELE_OFFSET_LIMIT)));

        for repeats in [6, 20] {
            let full = allele_support(RepeatCount(repeats));
            assert_eq!(
                full.len(),
                13,
                "at {repeats} repeats the support is the full ±{ALLELE_OFFSET_LIMIT}"
            );
            assert_eq!(full.first(), Some(&WholeRepeatOffset(-ALLELE_OFFSET_LIMIT)));
        }
    }

    /// A tract at the copy floors this caller routes on never sees this, but the arithmetic
    /// that clips the support must not panic on a repeat count that does not fit in the
    /// offset type — a satellite of 20,000 copies is one `u32` value away from any tract.
    #[test]
    fn the_support_survives_a_repeat_count_far_larger_than_an_offset() {
        let huge = allele_support(RepeatCount(u32::MAX));
        assert_eq!(huge.len(), 13);
        assert_eq!(huge.first(), Some(&WholeRepeatOffset(-ALLELE_OFFSET_LIMIT)));
    }

    /// A zero-copy tract is not something region typing emits, but the support must still be
    /// a set the fit can walk: at zero the only allele length that exists is the reference's
    /// own, and every longer one.
    #[test]
    fn a_tract_of_no_repeats_can_only_hold_alleles_at_or_above_the_reference() {
        let none = allele_support(RepeatCount(0));
        assert_eq!(none.first(), Some(&WholeRepeatOffset(0)));
        assert_eq!(
            none.len(),
            usize::try_from(ALLELE_OFFSET_LIMIT).unwrap() + 1
        );
    }
}
