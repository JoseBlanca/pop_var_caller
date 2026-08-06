//! What can go wrong with a read on the SNP/indel path, and what one cell is worth
//! under each genotype.
//!
//! **One thing can go wrong and nothing else: a base is misread.** A read over a
//! reference copy shows some other base with probability `ε`; a read over an alternative
//! copy reverts to the reference with probability `ε/3` — three bases to go wrong into,
//! one to come back to (`spec/parameter_prepass.md` §3). A repeat tract can also slip a
//! whole copy, which this model has nowhere to put, and that is why the STR path has its
//! own.
//!
//! **`ε` is broader than sequencing error**, and one sample cannot narrow it: it is the
//! rate at which a read shows an allele the individual does not carry, whatever put it
//! there — a misread base, a polymerase slip during amplification, a chimera, a damaged
//! base, contamination, index hopping, a read mismapped from a paralogue. The fitted
//! value is their sum, and the read-group grain fits the sum because most of them are
//! properties of a library and a run (`spec/parameter_prepass_generic.md` §2).
//!
//! **The part that is not obvious is the multi-library arm, and it is the reason this
//! file is its own commit.** A cell of the windowed table remembers how many alternative
//! reads came from each library and how many reads in total showed the reference — it
//! has forgotten how the *depth* split between libraries. The tempting repair is to
//! invent that split, giving each library its average share `n̂_g = w_g·n`. That is not a
//! probability over the cell space: it does not sum to one, and it charges a library for
//! reference reads it did not have whenever `k_g` exceeds its average share. Measured, it
//! reports heterozygosity **68% high** and the homozygous-non-reference rate **78% low**
//! at three reads a site — on two libraries with the *same* error rate, so the damage has
//! nothing to do with chemistry — and it does not shrink as data accumulates
//! (`doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md` §2). The rule
//! below sums over the forgotten split instead of guessing it, in closed form, at the
//! same cost.
//!
//! Design: `doc/devel/ng/arch/parameter_prepass_generic.md` §5.1,
//! `doc/devel/ng/spec/parameter_prepass_generic.md` §1 and §2.

use smallvec::SmallVec;

use crate::genetics::lgamma;
use crate::ng::parameter_estimation::fitting::NoiseModel;
use crate::ng::parameter_estimation::generic::histogram::{Attribution, Cell};
use crate::ng::types::{ErrorRate, Ploidy, ReadGroupId};

/// How far the shares of a sample's libraries may miss one before the set is refused.
///
/// They are computed as a division of read counts, so exact equality is not available;
/// this is loose enough for a sum of a handful of quotients and tight enough that a set
/// built from the wrong denominator cannot pass.
const SHARE_SUM_TOLERANCE: f64 = 1e-9;

/// How far below zero a cell's count of reference reads may sit before it is a fault
/// rather than a rounding.
///
/// A cell is scored at the mean of the exact depths that landed in it, and that mean is
/// never below the cell's own alternative count — `spec/parameter_prepass_generic.md`
/// §12.10 asserts it directly, and it is the one thing standing between this model and a
/// fit 5.2 rungs low with nothing on the outside to show for it. The two are equal when
/// every site in the cell had every read showing the alternative allele, and the division
/// that produced the mean can land a few ulp under.
const REFERENCE_READ_TOLERANCE: f64 = 1e-9;

/// One library's share of a sample's reads and how often one of its reads shows an
/// allele the individual does not carry.
///
/// The two numbers are named in a struct rather than travelling as a pair of `f64`s,
/// because a share and a rate are both small fractions and transposing them is a
/// plausible wrong answer rather than a compile error.
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct LibraryNoise {
    /// Which library.
    pub read_group: ReadGroupId,
    /// `w_g` — the share of this **sample's** reads this library produced, so the shares
    /// of a sample sum to one. A share, not a depth: the cell has forgotten the depths.
    pub share_of_reads: f64,
    /// `ε_g` — how often one of this library's reads shows an allele the individual does
    /// not carry.
    pub error_rate: ErrorRate,
}

/// Every library of one sample, with its share and its error rate: the SNP/indel path's
/// noise parameters, and what the profile scan steps through.
///
/// **The shares are in here beside the rates rather than in the model**, because the two
/// are inseparable in the expression and because parallel collections indexed by position
/// are how a library's rate gets scored against another library's share. The model itself
/// holds nothing.
///
/// Ascending by read group and with no duplicates, so a listing on a cell key — which is
/// canonical for the same reason — can be walked against it.
#[derive(Clone, PartialEq, Debug)]
pub struct SampleLibraryNoise {
    /// **Two inline.** 1,550 of the 1,707 samples in the tomato archive survey carry one
    /// library and the rest carry a handful, so the heap allocation is the exception
    /// (`spec/parameter_prepass_generic.md` §1).
    libraries: SmallVec<[LibraryNoise; 2]>,
}

impl SampleLibraryNoise {
    /// The libraries of one sample, in any order.
    ///
    /// # Panics
    ///
    /// If the sample has no libraries, if a read group appears twice, if a share is not a
    /// fraction strictly above zero, or if the shares do not sum to one. Each of those is
    /// a caller that built the set from the wrong reads, and each would otherwise produce
    /// a score that is not a probability over the cell space — which is precisely the
    /// failure this file exists to avoid.
    #[must_use]
    pub fn new(libraries: impl IntoIterator<Item = LibraryNoise>) -> Self {
        let mut libraries: SmallVec<[LibraryNoise; 2]> = libraries.into_iter().collect();
        assert!(
            !libraries.is_empty(),
            "a sample with no libraries has no reads to have an error rate"
        );
        libraries.sort_unstable_by_key(|library| library.read_group);

        let mut total_share = 0.0;
        for (position, library) in libraries.iter().enumerate() {
            assert!(
                library.share_of_reads > 0.0 && library.share_of_reads <= 1.0,
                "read group {} produced share {} of the sample's reads, which is not a \
                 fraction of them",
                library.read_group.get(),
                library.share_of_reads
            );
            assert!(
                position == 0 || libraries[position - 1].read_group != library.read_group,
                "read group {} is listed twice",
                library.read_group.get()
            );
            total_share += library.share_of_reads;
        }
        assert!(
            (total_share - 1.0).abs() < SHARE_SUM_TOLERANCE,
            "the libraries' shares sum to {total_share} rather than one"
        );

        Self { libraries }
    }

    /// A sample sequenced from one library, which is the common case and the whole of the
    /// read-group table: every entry there belongs to one group already, so it is scored
    /// as though that group were the sample.
    #[must_use]
    pub fn single(read_group: ReadGroupId, error_rate: ErrorRate) -> Self {
        Self {
            libraries: SmallVec::from_elem(
                LibraryNoise {
                    read_group,
                    share_of_reads: 1.0,
                    error_rate,
                },
                1,
            ),
        }
    }

    /// The libraries, ascending by read group.
    #[must_use]
    pub fn libraries(&self) -> &[LibraryNoise] {
        &self.libraries
    }

    /// One library by its read group.
    ///
    /// # Panics
    ///
    /// If the sample has no such library. A cell listing a read group these parameters do
    /// not know is a cell and a parameter set built from two different read-group tables,
    /// and scoring it against some other library's rate is exactly the mistake the
    /// attributed key exists to prevent.
    fn library(&self, read_group: ReadGroupId) -> &LibraryNoise {
        let position = self
            .libraries
            .binary_search_by_key(&read_group, |library| library.read_group)
            .unwrap_or_else(|_| {
                panic!(
                    "a cell lists read group {} and these noise parameters do not have it",
                    read_group.get()
                )
            });
        &self.libraries[position]
    }
}

/// The SNP/indel path's noise model: a per-base substitution rate and nothing else.
///
/// Stateless — everything it needs travels in [`SampleLibraryNoise`] — so that one value
/// serves every scan.
#[derive(Copy, Clone, PartialEq, Eq, Debug, Default)]
pub struct SubstitutionNoiseModel;

impl NoiseModel for SubstitutionNoiseModel {
    type Cell = Cell;
    type NoiseParams = SampleLibraryNoise;

    /// `ln L(cell | genotype)` for every genotype, appended in ascending order of
    /// alternative copies.
    ///
    /// Each read independently picks a library — library `g` with probability `w_g` —
    /// and then shows the alternative allele with probability `p_j(ε_g)` or the reference
    /// with `1 − p_j(ε_g)`. So a cell of the **attributed** arm, which records how many
    /// alternative reads came from each library and how many reads showed the reference
    /// in total, is one multinomial over `G + 1` categories:
    ///
    /// ```text
    ///                        n!                                                   n−k
    /// L(cell | j)  =  ────────────────  ·  Π (w_g·p_j(ε_g))^{k_g}  ·  ( Σ w_g·(1 − p_j(ε_g)) )
    ///                  Π k_g! (n−k)!       g                           g
    /// ```
    ///
    /// **The sum in the last factor runs over every library the sample has**, listed at
    /// this cell or not — a library that showed no alternative read here still produced
    /// reads that showed the reference. The product before it runs only over the ones
    /// listed, which is the same thing, since a `k_g` of zero contributes a factor of one.
    ///
    /// A **pooled** cell is the same expression with the `G` alternative categories
    /// collapsed into one, which leaves a binomial at the share-weighted rate
    /// `Σ_g w_g·p_j(ε_g)`. Nothing is approximated in either: the pooled form is the exact
    /// sum of the attributed form over the splits the pooled key forgot.
    ///
    /// **The depth is a real number**, because a binned cell is scored at the mean of the
    /// exact depths that fell in it. Everything but the two factorial prefactors is
    /// already continuous in the depth, and those two are common to every genotype, so
    /// they cancel out of the mixture and cannot move a fit — they are computed exactly
    /// anyway, because the identity that the rule sums to one over the cell space is what
    /// separates it from the plug-in it replaced.
    ///
    /// # Panics
    ///
    /// If `ploidy` disagrees with the cell's own, if the cell is charged a negative count
    /// of reference reads, or if it lists a read group these parameters do not have.
    fn append_genotype_likelihoods(
        &self,
        cell: &Cell,
        noise: &SampleLibraryNoise,
        ploidy: Ploidy,
        out: &mut Vec<f64>,
    ) {
        assert_eq!(
            ploidy, cell.ploidy,
            "a cell of ploidy {:?} scored against the genotypes of ploidy {ploidy:?}",
            cell.ploidy
        );

        let depth = cell.mean_depth;
        let alt_reads = cell.key.alt_reads();
        let reference_reads = depth - f64::from(alt_reads);
        // `spec/parameter_prepass_generic.md` §12.8's second identity, asserted where it
        // would first do damage rather than left to the accumulator. A negative exponent
        // here is a cell scored at a depth below its own alternative count, which is what
        // taking the mean over the whole bin instead of the cell produces on 0.3% of
        // sites — and it lands the fit 5.2 rungs low with nothing to show for it.
        assert!(
            reference_reads >= -REFERENCE_READ_TOLERANCE,
            "a cell of {alt_reads} alternative reads is scored at depth {depth}, which \
             charges it {reference_reads} reference reads"
        );
        let reference_reads = reference_reads.max(0.0);

        out.reserve(usize::from(ploidy.get()) + 1);
        for alt_copies in 0..=ploidy.get() {
            out.push(match cell.key.attribution() {
                Attribution::Pooled => ln_likelihood_pooled(
                    depth,
                    alt_reads,
                    reference_reads,
                    noise,
                    ploidy,
                    alt_copies,
                ),
                Attribution::ByReadGroup(listing) => ln_likelihood_attributed(
                    depth,
                    reference_reads,
                    listing,
                    noise,
                    ploidy,
                    alt_copies,
                ),
            });
        }
    }
}

/// `p_j(ε)` — the chance **one read** shows something other than the reference base.
///
/// It draws one of the individual's `P` copies, hits an alternative copy with probability
/// `j/P`, and either reading can be misread: an alternative copy still reads
/// non-reference unless it is misread into the one reference base, and a reference copy
/// reads non-reference if it is misread into any of the three others
/// (`spec/parameter_prepass.md` §3).
///
/// At `P = 2` and `j = 1` this is `½ + ε/3`, not `½`: a misread of the reference copy
/// *adds* a non-reference read at `ε/2` while a misread of the alternative copy *removes*
/// one at only `ε/6`, so the expected non-reference fraction at a heterozygote sits above
/// a half. Three parts in ten thousand at `ε = 0.001`, and free to get right.
fn alternative_read_probability(alt_copies: u8, ploidy: Ploidy, error_rate: f64) -> f64 {
    let carried = f64::from(alt_copies) / f64::from(ploidy.get());
    carried * (1.0 - error_rate / 3.0) + (1.0 - carried) * error_rate
}

/// The share-weighted rate `Σ_g w_g·p_j(ε_g)` a pooled cell's alternative reads are drawn
/// at, and the rate `Σ_g w_g·(1 − p_j(ε_g))` its reference reads are.
///
/// **The two are exactly complementary, and only because the shares sum to one** —
/// `Σ_g w_g·(1 − p_g)` is `Σ_g w_g − Σ_g w_g·p_g`, which is `1 − Σ_g w_g·p_g`. Writing
/// either in place of the other changes no answer, which a mutation confirmed. Both are
/// summed out longhand anyway, because that is the form
/// `spec/parameter_prepass_generic.md` §1 states and the shares' sum is an invariant of
/// [`SampleLibraryNoise`] rather than of this expression.
///
/// What is *not* interchangeable is which libraries each runs over. The reference rate is
/// summed over **every** library the sample has, listed at the cell or not: a library that
/// showed no alternative read here still produced reads that showed the reference. That
/// is the term the retired plug-in got wrong.
fn share_weighted_rates(noise: &SampleLibraryNoise, ploidy: Ploidy, alt_copies: u8) -> (f64, f64) {
    let mut alternative = 0.0;
    let mut reference = 0.0;
    for library in noise.libraries() {
        let p = alternative_read_probability(alt_copies, ploidy, library.error_rate.get());
        alternative += library.share_of_reads * p;
        reference += library.share_of_reads * (1.0 - p);
    }
    (alternative, reference)
}

/// A cell whose key did not keep which library the alternative reads came from: one
/// binomial at the share-weighted rate.
fn ln_likelihood_pooled(
    depth: f64,
    alt_reads: u32,
    reference_reads: f64,
    noise: &SampleLibraryNoise,
    ploidy: Ploidy,
    alt_copies: u8,
) -> f64 {
    let (alternative_rate, reference_rate) = share_weighted_rates(noise, ploidy, alt_copies);
    ln_factorial_real(depth) - ln_factorial(alt_reads) - ln_factorial_real(reference_reads)
        + count_times_ln(f64::from(alt_reads), alternative_rate)
        + count_times_ln(reference_reads, reference_rate)
}

/// A cell whose key kept which library each alternative read came from: one multinomial
/// over one category per library plus a pooled "showed the reference".
fn ln_likelihood_attributed(
    depth: f64,
    reference_reads: f64,
    listing: &[(ReadGroupId, u8)],
    noise: &SampleLibraryNoise,
    ploidy: Ploidy,
    alt_copies: u8,
) -> f64 {
    let (_, reference_rate) = share_weighted_rates(noise, ploidy, alt_copies);
    let mut total = ln_factorial_real(depth) - ln_factorial_real(reference_reads);
    for &(read_group, alt_from_group) in listing {
        let library = noise.library(read_group);
        let p = alternative_read_probability(alt_copies, ploidy, library.error_rate.get());
        total += -ln_factorial(u32::from(alt_from_group))
            + count_times_ln(f64::from(alt_from_group), library.share_of_reads * p);
    }
    total + count_times_ln(reference_reads, reference_rate)
}

/// `count · ln probability`, with the convention that a count of zero contributes
/// nothing.
///
/// Without it a category no read fell into, at a probability of zero, would be
/// `0 · −∞ = NaN` — and a `NaN` does not lose loudly in a scan, it simply never wins.
fn count_times_ln(count: f64, probability: f64) -> f64 {
    if count == 0.0 {
        0.0
    } else {
        count * probability.ln()
    }
}

/// `ln n!` for a whole count.
fn ln_factorial(n: u32) -> f64 {
    lgamma(f64::from(n) + 1.0)
}

/// `ln x!` for a **fractional** `x`, which is what a cell scored at a mean depth needs.
///
/// Defined through `ln Γ(x + 1)`, so `x` must be above `−1`. A cell cannot reach that:
/// its depth is a mean of positive depths and its reference-read count is guarded above.
fn ln_factorial_real(x: f64) -> f64 {
    lgamma(x + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::generic::depth_bins::DepthBinEdges;
    use crate::ng::parameter_estimation::generic::histogram::SiteKey;

    fn group(id: u32) -> ReadGroupId {
        ReadGroupId(id)
    }

    fn rate(value: f64) -> ErrorRate {
        ErrorRate::try_new(value).expect("a probability")
    }

    fn ploidy(copies: u8) -> Ploidy {
        Ploidy::try_new(copies).expect("a positive copy number")
    }

    /// Two libraries with a four-fold difference in error rate and a 90/10 split of the
    /// reads — one of the research harness's own worlds, and the regime the attributed
    /// key exists for.
    fn two_libraries() -> SampleLibraryNoise {
        SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 0.9,
                error_rate: rate(1e-3),
            },
            LibraryNoise {
                read_group: group(1),
                share_of_reads: 0.1,
                error_rate: rate(4e-3),
            },
        ])
    }

    /// A cell of the given depth and library breakdown, at the given ploidy. The bin is
    /// the one the depth falls in, which the scoring rule never reads — it reads
    /// `mean_depth` — but a key needs one and a wrong one would misfile the cell.
    fn cell_of(depth: f64, alt_by_group: &[(ReadGroupId, u32)], copies: u8) -> Cell {
        let edges = DepthBinEdges::new();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let bin = edges.bin_for(depth.round() as u32);
        Cell {
            key: SiteKey::attributing(bin, alt_by_group),
            ploidy: ploidy(copies),
            sites: 1,
            mean_depth: depth,
        }
    }

    /// A pooled cell — every entry of the read-group table, and every entry of a
    /// single-library sample's windowed one.
    fn pooled_cell_of(depth: f64, alt_reads: u32, copies: u8) -> Cell {
        let edges = DepthBinEdges::new();
        #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
        let bin = edges.bin_for(depth.round() as u32);
        Cell {
            key: SiteKey::pooled(bin, alt_reads),
            ploidy: ploidy(copies),
            sites: 1,
            mean_depth: depth,
        }
    }

    fn score(cell: &Cell, noise: &SampleLibraryNoise) -> Vec<f64> {
        let mut out = Vec::new();
        SubstitutionNoiseModel.append_genotype_likelihoods(cell, noise, cell.ploidy, &mut out);
        out
    }

    /// Every cell a sample of `libraries` libraries at depth `depth` can produce, keyed
    /// the way the accumulator keys them: the library breakdown is kept while the total
    /// alternative count is at most [`MAX_ATTRIBUTED_ALT_READS`], and pooled above.
    ///
    /// [`MAX_ATTRIBUTED_ALT_READS`]: crate::ng::parameter_estimation::generic::histogram::MAX_ATTRIBUTED_ALT_READS
    fn whole_cell_space(depth: u32, libraries: usize, copies: u8) -> Vec<Cell> {
        use crate::ng::parameter_estimation::generic::histogram::MAX_ATTRIBUTED_ALT_READS;

        let mut cells = Vec::new();
        for alt_reads in 0..=depth {
            if alt_reads == 0 || alt_reads > MAX_ATTRIBUTED_ALT_READS {
                cells.push(pooled_cell_of(f64::from(depth), alt_reads, copies));
                continue;
            }
            for split in compositions(alt_reads, libraries) {
                let listing: Vec<(ReadGroupId, u32)> = split
                    .iter()
                    .enumerate()
                    .map(|(library, &count)| (group(library as u32), count))
                    .collect();
                cells.push(cell_of(f64::from(depth), &listing, copies));
            }
        }
        cells
    }

    /// Every way of splitting `total` among `parts` libraries, zeros included.
    fn compositions(total: u32, parts: usize) -> Vec<Vec<u32>> {
        if parts == 1 {
            return vec![vec![total]];
        }
        let mut out = Vec::new();
        for first in 0..=total {
            for rest in compositions(total - first, parts - 1) {
                let mut one = vec![first];
                one.extend(rest);
                out.push(one);
            }
        }
        out
    }

    /// **Identity 1 of `spec/parameter_prepass_generic.md` §12.8: the rule is a
    /// probability over the cell space, at any parameter values.**
    ///
    /// Every cell a site of this depth can land in, keyed exactly as the accumulator
    /// keys it — attributed below the bound and pooled above — summed under one genotype,
    /// must come to one. It is the check the average-share plug-in fails, and it fails it
    /// the day it is written: a rule that is not a probability cannot be unbiased.
    ///
    /// Both arms are exercised at every depth here, since 5 is above the bound of four.
    #[test]
    fn the_rule_sums_to_one_over_the_cell_space() {
        let noise = two_libraries();
        for depth in [1_u32, 3, 5, 9] {
            for copies in [2_u8, 4] {
                let cells = whole_cell_space(depth, 2, copies);
                for alt_copies in 0..=copies {
                    let total: f64 = cells
                        .iter()
                        .map(|cell| score(cell, &noise)[usize::from(alt_copies)].exp())
                        .sum();
                    assert!(
                        (total - 1.0).abs() < 1e-12,
                        "depth {depth}, ploidy {copies}, {alt_copies} alternative \
                         copies: the cell space sums to {total}"
                    );
                }
            }
        }
    }

    /// The same identity on four libraries, where the attributed arm has 35 splits of
    /// four alternative reads rather than five and the multinomial is doing real work.
    #[test]
    fn the_rule_sums_to_one_over_a_four_library_cell_space() {
        let noise = SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 0.25,
                error_rate: rate(5e-4),
            },
            LibraryNoise {
                read_group: group(1),
                share_of_reads: 0.25,
                error_rate: rate(1e-3),
            },
            LibraryNoise {
                read_group: group(2),
                share_of_reads: 0.25,
                error_rate: rate(2e-3),
            },
            LibraryNoise {
                read_group: group(3),
                share_of_reads: 0.25,
                error_rate: rate(4e-3),
            },
        ]);

        for depth in [2_u32, 6] {
            let cells = whole_cell_space(depth, 4, 2);
            for alt_copies in 0..=2_usize {
                let total: f64 = cells
                    .iter()
                    .map(|cell| score(cell, &noise)[alt_copies].exp())
                    .sum();
                assert!(
                    (total - 1.0).abs() < 1e-12,
                    "depth {depth}, {alt_copies} alternative copies: sums to {total}"
                );
            }
        }
    }

    /// The identity has to hold **at any parameters**, not only at plausible ones — it is
    /// an algebraic property of the expression and not a numerical accident near the
    /// error rates a sequencer produces. A rate of 0.3 is far above anything the ladder
    /// reaches.
    #[test]
    fn the_rule_sums_to_one_at_implausible_error_rates_too() {
        let noise = SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 0.5,
                error_rate: rate(0.3),
            },
            LibraryNoise {
                read_group: group(1),
                share_of_reads: 0.5,
                error_rate: rate(1e-5),
            },
        ]);
        let cells = whole_cell_space(6, 2, 2);
        for alt_copies in 0..=2_usize {
            let total: f64 = cells
                .iter()
                .map(|cell| score(cell, &noise)[alt_copies].exp())
                .sum();
            assert!((total - 1.0).abs() < 1e-12, "sums to {total}");
        }
    }

    /// **Identity 2: no cell is charged a negative count of reference reads.** The bound
    /// this rests on is `mean_depth ≥ alt_reads`, which holds by construction when a cell
    /// is scored at the mean of its **own** depths and fails on 0.3% of sites when the
    /// mean is taken over the whole bin.
    ///
    /// Asserted where it does damage: over the deep, alternative-rich corner of the
    /// table, which is the exact case the per-bin mean fails.
    #[test]
    fn no_cell_is_charged_a_negative_count_of_reference_reads() {
        let noise = two_libraries();
        for depth in 100..=124_u32 {
            for alt_reads in 0..=depth {
                let cell = pooled_cell_of(f64::from(depth), alt_reads, 2);
                let scores = score(&cell, &noise);
                assert_eq!(scores.len(), 3);
                assert!(
                    scores.iter().all(|value| *value <= 0.0),
                    "depth {depth}, {alt_reads} alternative reads: {scores:?}"
                );
            }
        }
    }

    /// A cell whose depth sits below its own alternative count is not scored, it is
    /// named — because the number it would produce is plausible and bounded, so nothing
    /// downstream would show it.
    #[test]
    #[should_panic(expected = "charges it -3 reference reads")]
    fn a_cell_scored_below_its_own_alternative_count_is_refused() {
        let noise = two_libraries();
        let cell = pooled_cell_of(9.0, 12, 2);
        let _ = score(&cell, &noise);
    }

    /// **Identity 3: with every library's error rate equal there is nothing to
    /// attribute, so the attributed rule must reproduce the pooled one exactly.**
    ///
    /// Not "up to a constant": at one library the two expressions are the same
    /// expression, since `w = 1` collapses the multinomial coefficient back to a
    /// binomial one. That is the sharper half of the check and it is asserted first.
    #[test]
    fn at_one_library_the_attributed_rule_is_the_pooled_rule() {
        let noise = SampleLibraryNoise::single(group(0), rate(1e-3));
        for depth in [5.0_f64, 7.5, 40.0] {
            for alt_reads in 1..=4_u32 {
                let attributed = score(&cell_of(depth, &[(group(0), alt_reads)], 2), &noise);
                let pooled = score(&pooled_cell_of(depth, alt_reads, 2), &noise);
                for (alt_copies, (&got, &want)) in attributed.iter().zip(&pooled).enumerate() {
                    assert!(
                        (got - want).abs() < 1e-12,
                        "depth {depth}, {alt_reads} alternative reads, {alt_copies} \
                         copies: attributed {got}, pooled {want}"
                    );
                }
            }
        }
    }

    /// The multi-library half of identity 3. With the rates equal the attributed rule
    /// still keeps a per-cell factor the pooled one does not — the probability of *that*
    /// split of the alternative reads — but it carries no genotype, so it cancels out of
    /// the mixture and the estimator cannot see it. What must agree exactly, then, is the
    /// **differences between genotypes**, which is what a fit is computed from.
    #[test]
    fn at_equal_error_rates_the_attributed_and_pooled_rules_differ_by_no_genotype() {
        let flat = SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 0.9,
                error_rate: rate(2e-3),
            },
            LibraryNoise {
                read_group: group(1),
                share_of_reads: 0.1,
                error_rate: rate(2e-3),
            },
        ]);

        for depth in [4.0_f64, 11.0, 60.0] {
            for split in [
                vec![(group(0), 1_u32)],
                vec![(group(1), 1)],
                vec![(group(0), 2), (group(1), 1)],
                vec![(group(0), 1), (group(1), 3)],
            ] {
                let alt_reads: u32 = split.iter().map(|&(_, count)| count).sum();
                let attributed = score(&cell_of(depth, &split, 2), &flat);
                let pooled = score(&pooled_cell_of(depth, alt_reads, 2), &flat);
                for alt_copies in 1..=2_usize {
                    let got = attributed[alt_copies] - attributed[0];
                    let want = pooled[alt_copies] - pooled[0];
                    assert!(
                        (got - want).abs() < 1e-12,
                        "depth {depth}, split {split:?}, {alt_copies} copies: \
                         attributed {got}, pooled {want}"
                    );
                }
            }
        }
    }

    /// And the same difference must **not** vanish when the rates differ, or the
    /// attributed key would be keeping something it cannot use — which is the whole claim
    /// the key rests on. The pooled key sees only the share-weighted mean `ε̄`, so on a
    /// four-fold ratio the two arms have to disagree.
    #[test]
    fn at_unequal_error_rates_which_library_the_reads_came_from_changes_the_score() {
        let noise = two_libraries();
        let cheap = score(&cell_of(6.0, &[(group(0), 2)], 2), &noise);
        let dear = score(&cell_of(6.0, &[(group(1), 2)], 2), &noise);

        let from_cheap = cheap[1] - cheap[0];
        let from_dear = dear[1] - dear[0];
        assert!(
            (from_cheap - from_dear).abs() > 0.1,
            "two alternative reads from the 0.001 library scored {from_cheap} against \
             {from_dear} from the 0.004 one; the key is keeping nothing"
        );
    }

    /// The pooled arm is the attributed arm summed over the splits the pooled key forgot,
    /// exactly — which is why the two can sit in one table without the total ceasing to
    /// be a probability. Checked directly rather than inferred from the sum-to-one test.
    #[test]
    fn a_pooled_cell_is_its_attributed_cells_added_up() {
        let noise = two_libraries();
        for alt_reads in 1..=4_u32 {
            let pooled = score(&pooled_cell_of(8.0, alt_reads, 2), &noise);
            for (alt_copies, &pooled_score) in pooled.iter().enumerate() {
                let summed: f64 = compositions(alt_reads, 2)
                    .into_iter()
                    .map(|split| {
                        let listing: Vec<(ReadGroupId, u32)> = split
                            .iter()
                            .enumerate()
                            .map(|(library, &count)| (group(library as u32), count))
                            .collect();
                        score(&cell_of(8.0, &listing, 2), &noise)[alt_copies].exp()
                    })
                    .sum();
                assert!(
                    (summed - pooled_score.exp()).abs() < 1e-15,
                    "{alt_reads} alternative reads, {alt_copies} copies: splits sum to \
                     {summed}, pooled is {}",
                    pooled_score.exp()
                );
            }
        }
    }

    /// **Identity 4: the same numbers the research harness computes.**
    ///
    /// The rows below are `ln_component_attributed`'s own output on the harness's
    /// `ratio=4 depth=6 split=skew90` world — two libraries at `ε` = 0.001 and 0.004 with
    /// a 90/10 split — generated by
    /// `cargo run --release --example ng_multilib_key_harness -- --only=oracle`. That
    /// program is where the multi-library key was measured to be unbiased in all 31
    /// worlds tried, so agreeing with it is agreeing with the thing the design was
    /// decided on
    /// (`doc/devel/ng/research/parameter_estimator_experiments_2026-08-06.md` §2).
    ///
    /// **Compared as differences between genotypes, and that is the sharper comparison
    /// rather than the looser one.** The two implementations reach `ln Γ` differently —
    /// the harness by a Lanczos series, this file through `libm`'s — so their factorial
    /// prefactors disagree in the last bits. Those prefactors carry no genotype, so they
    /// cancel from every difference, and what is left is exactly what a fit can see: two
    /// rules that agree on the differences give the same responsibilities, the same
    /// climb, and the same answer. The absolute value of the prefactors is held instead
    /// by `the_rule_sums_to_one_over_the_cell_space`, which cannot pass if they are wrong.
    #[test]
    fn the_rule_agrees_with_the_research_harness_on_one_of_its_worlds() {
        const HARNESS: [(f64, [u32; 2], [f64; 3]); 18] = [
            (
                3.0,
                [1, 0],
                [
                    -5.91710519743795e0,
                    -1.0872574090050942e0,
                    -1.4495088222256634e1,
                ],
            ),
            (
                3.0,
                [0, 1],
                [
                    -6.728035413654278e0,
                    -3.2824853134649903e0,
                    -1.6693313633704754e1,
                ],
            ),
            (
                3.0,
                [2, 1],
                [
                    -2.075166531146811e1,
                    -3.49013937081562e0,
                    -1.4166948364572005e0,
                ],
            ),
            (
                6.0,
                [1, 0],
                [
                    -5.227860554077148e0,
                    -2.476152897443037e0,
                    -3.70339509512452e1,
                ],
            ),
            (
                6.0,
                [0, 1],
                [
                    -6.038790770293477e0,
                    -4.6713808019029335e0,
                    -3.923217636269332e1,
                ],
            ),
            (
                6.0,
                [2, 1],
                [
                    -1.7759835575113264e1,
                    -2.576449766259516e0,
                    -2.165297247245173e1,
                ],
            ),
            (
                6.0,
                [1, 3],
                [
                    -2.6393510956452833e1,
                    -7.0707326038546245e0,
                    -1.841111389672419e1,
                ],
            ),
            (
                6.0,
                [0, 4],
                [
                    -2.859073553378905e1,
                    -1.0652254869434408e1,
                    -2.1995633669292197e1,
                ],
            ),
            (
                6.5,
                [1, 0],
                [
                    -5.14846826927014e0,
                    -2.7431173012691525e0,
                    -4.08259098951631e1,
                ],
            ),
            (
                6.5,
                [0, 1],
                [
                    -5.959398485486468e0,
                    -4.938345205729048e0,
                    -4.3024135306611214e1,
                ],
            ),
            (
                6.5,
                [2, 1],
                [
                    -1.7467350074845545e1,
                    -2.630320954624922e0,
                    -2.5231838200908907e1,
                ],
            ),
            (
                6.5,
                [1, 3],
                [
                    -2.5946874776357856e1,
                    -6.9704531123927715e0,
                    -2.1835828945354105e1,
                ],
            ),
            (
                6.5,
                [0, 4],
                [
                    -2.8144099353694077e1,
                    -1.0551975377972555e1,
                    -2.542034871792211e1,
                ],
            ),
            (
                20.0,
                [1, 0],
                [
                    -4.042099590013885e0,
                    -1.098837921510725e1,
                    -1.44246024391479e2,
                ],
            ),
            (
                20.0,
                [0, 1],
                [
                    -4.8530298062302135e0,
                    -1.3183607119567146e1,
                    -1.4644424980292715e2,
                ],
            ),
            (
                20.0,
                [2, 1],
                [
                    -1.3734996147541397e1,
                    -8.249597620415127e0,
                    -1.2602596744917695e2,
                ],
            ),
            (
                20.0,
                [1, 3],
                [
                    -2.0634070473492855e1,
                    -1.1009279402622123e1,
                    -1.210495078180613e2,
                ],
            ),
            (
                20.0,
                [0, 4],
                [
                    -2.2831295050829077e1,
                    -1.4590801668201907e1,
                    -1.246340275906293e2,
                ],
            ),
        ];

        let noise = two_libraries();
        for (depth, split, want) in HARNESS {
            let listing: Vec<(ReadGroupId, u32)> = vec![(group(0), split[0]), (group(1), split[1])];
            let got = score(&cell_of(depth, &listing, 2), &noise);
            for alt_copies in 1..=2_usize {
                let ours = got[alt_copies] - got[0];
                let theirs = want[alt_copies] - want[0];
                assert!(
                    (ours - theirs).abs() < 1e-12,
                    "depth {depth}, split {split:?}, {alt_copies} copies: ours {ours}, \
                     the harness's {theirs}"
                );
            }
        }
    }

    /// `p_j` at the three diploid rows of `spec/parameter_prepass.md` §3's table, which
    /// is where the `3` lives: a reference copy goes wrong into any of three other bases
    /// at `ε`, an alternative copy comes back into the one reference base at `ε/3`.
    /// Writing `1 − ε` in the last row charges the alternative copy three times too much
    /// reversion and biases the homozygous-non-reference rate.
    #[test]
    fn the_per_read_probability_is_the_specs_three_rows() {
        let error_rate = 0.003;
        let diploid = ploidy(2);
        assert!((alternative_read_probability(0, diploid, error_rate) - error_rate).abs() < 1e-15);
        assert!(
            (alternative_read_probability(1, diploid, error_rate) - (0.5 + error_rate / 3.0)).abs()
                < 1e-15
        );
        assert!(
            (alternative_read_probability(2, diploid, error_rate) - (1.0 - error_rate / 3.0)).abs()
                < 1e-15
        );

        // A tetraploid's intermediate dosages are the same line at j/P.
        let tetraploid = ploidy(4);
        let quarter = alternative_read_probability(1, tetraploid, error_rate);
        assert!((quarter - (0.25 * (1.0 - error_rate / 3.0) + 0.75 * error_rate)).abs() < 1e-15);
    }

    /// The appending contract the profile scan is built on: two cells in a row leave one
    /// buffer holding both rows, in order, with nothing cleared.
    #[test]
    fn the_model_appends_rather_than_clearing() {
        let noise = two_libraries();
        let mut out = vec![f64::NAN];
        let first = pooled_cell_of(5.0, 0, 2);
        let second = pooled_cell_of(5.0, 1, 2);

        SubstitutionNoiseModel.append_genotype_likelihoods(&first, &noise, ploidy(2), &mut out);
        SubstitutionNoiseModel.append_genotype_likelihoods(&second, &noise, ploidy(2), &mut out);

        assert_eq!(out.len(), 7, "one sentinel and two rows of three");
        assert!(out[0].is_nan(), "the buffer's prior contents were cleared");
        assert_eq!(&out[1..4], score(&first, &noise).as_slice());
        assert_eq!(&out[4..7], score(&second, &noise).as_slice());
    }

    /// A tetraploid cell gets five genotypes, not three — the loop bound is the cell's
    /// own ploidy and nothing here is diploid.
    #[test]
    fn a_tetraploid_cell_is_scored_against_five_genotypes() {
        let noise = two_libraries();
        let scores = score(&pooled_cell_of(12.0, 3, 4), &noise);
        assert_eq!(scores.len(), 5);
        assert!(scores.iter().all(|value| value.is_finite()));
    }

    #[test]
    #[should_panic(
        expected = "a cell of ploidy Ploidy(4) scored against the genotypes of ploidy Ploidy(2)"
    )]
    fn a_cell_scored_against_another_ploidys_genotypes_is_refused() {
        let noise = two_libraries();
        let cell = pooled_cell_of(12.0, 3, 4);
        let mut out = Vec::new();
        SubstitutionNoiseModel.append_genotype_likelihoods(&cell, &noise, ploidy(2), &mut out);
    }

    #[test]
    #[should_panic(
        expected = "a cell lists read group 3 and these noise parameters do not have it"
    )]
    fn a_cell_naming_a_library_the_parameters_do_not_have_is_refused() {
        let noise = two_libraries();
        let cell = cell_of(6.0, &[(group(3), 2)], 2);
        let _ = score(&cell, &noise);
    }

    #[test]
    #[should_panic(expected = "the libraries' shares sum to 0.8 rather than one")]
    fn libraries_whose_shares_do_not_sum_to_one_are_refused() {
        let _ = SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 0.5,
                error_rate: rate(1e-3),
            },
            LibraryNoise {
                read_group: group(1),
                share_of_reads: 0.3,
                error_rate: rate(1e-3),
            },
        ]);
    }

    #[test]
    #[should_panic(expected = "read group 0 is listed twice")]
    fn a_library_listed_twice_is_refused() {
        let _ = SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 0.5,
                error_rate: rate(1e-3),
            },
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 0.5,
                error_rate: rate(4e-3),
            },
        ]);
    }

    #[test]
    #[should_panic(expected = "produced share 0 of the sample's reads")]
    fn a_library_that_produced_no_reads_is_refused() {
        let _ = SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(0),
                share_of_reads: 1.0,
                error_rate: rate(1e-3),
            },
            LibraryNoise {
                read_group: group(1),
                share_of_reads: 0.0,
                error_rate: rate(4e-3),
            },
        ]);
    }

    #[test]
    #[should_panic(expected = "a sample with no libraries")]
    fn a_sample_with_no_libraries_is_refused() {
        let _ = SampleLibraryNoise::new([]);
    }

    /// The listing is canonical however it arrives, which the cell keys rely on.
    #[test]
    fn the_libraries_come_out_ascending_by_read_group() {
        let noise = SampleLibraryNoise::new([
            LibraryNoise {
                read_group: group(7),
                share_of_reads: 0.4,
                error_rate: rate(1e-3),
            },
            LibraryNoise {
                read_group: group(2),
                share_of_reads: 0.6,
                error_rate: rate(4e-3),
            },
        ]);
        let ids: Vec<u32> = noise
            .libraries()
            .iter()
            .map(|library| library.read_group.get())
            .collect();
        assert_eq!(ids, vec![2, 7]);
    }
}
