//! **The denominator of the read likelihood's error-rate scale, per read group.**
//!
//! The caller charges each read the error probability the walk minted for it, rescaled by one
//! number per read group so that the average over that read group's admitted reads comes out at
//! the rate the pre-pass measured
//! (`doc/devel/ng/spec/read_likelihoods.md` §3.2):
//!
//! ```text
//!                fitted error rate for this read group
//! scale  =  ────────────────────────────────────────────────
//!           average minted error over that group's own reads
//! ```
//!
//! The pre-pass fits the numerator. This module carries the denominator, and nothing else.
//!
//! # It sums, and does not mint
//!
//! Every number here is one the walk has already produced. A read's minted error is the worse of
//! its window's base quality and its mapping quality, in log space
//! ([`minted_ln_read_error`](crate::ng::locus_generation::pileup)); the fold sums those into an
//! observation's `q_sum` and keeps the count in `num_obs`, and this adds those two up per read
//! group. **So there is no second definition of "how wrong is this read" to drift from the first**,
//! which is the requirement §3.2 spends most of its length on, and no traversal of the genome that
//! was not already being made.
//!
//! # The average is the geometric mean, and the spec used to ask for the other one
//!
//! `q_sum` is a sum of *logarithms*. So `exp(Σ q_sum / Σ num_obs)` is the geometric mean of the
//! minted error probabilities, and the arithmetic mean cannot be recovered from it — the walk
//! throws the individual reads away.
//!
//! **That is the right average anyway, and the reason is what the scale is applied to.** The model
//! charges an observation `exp(q_sum / num_obs)`, and so does production
//! (`var_calling/posterior_engine.rs`, which has no recalibration at all — there is nothing to copy
//! there but the quantity). A scale built from an arithmetic mean and applied to a geometric one
//! would not make the calibrated property hold in the model's own terms, so supplying the
//! arithmetic sum would have bought an inexactness rather than removed one — at the price of a
//! second accumulation at fold time and a field on every observation. Corrected in the spec on
//! 2026-08-24 by the owner.
//!
//! # Why the per-position depth cap does not divide the two site sets
//!
//! §3.2 requires the average to run over exactly the reads the fitted rate was fitted from, and the
//! histogram route thins every position to at most
//! [`MAX_BINNED_DEPTH`](super::depth_bins) reads before fitting. That would matter for a sum and
//! does not matter for a mean: the draw is hypergeometric on counts and never looks at a read's
//! quality, so the mean log error over the kept reads has the same expectation as over all of them.
//! **What has to match is the count against the mean it belongs to**, and both come from the same
//! observations.

use std::collections::BTreeMap;

use crate::ng::locus_generation::{LocusKind, SampleLocusObservations};
use crate::ng::types::ReadGroupId;

/// How many parts of one the log-error sum is counted in — see [`MintedReadErrors`].
const PARTS_OF_ONE: f64 = (1_u64 << 20) as f64;

/// One read group's minted error, summed over the sites a fit read.
///
/// **Two scalars and not a mean**, because a mean cannot be added to another shard's or another
/// sample's. The pre-pass accumulates in region shards and fits per sample, and the scale is
/// wanted per read group over a whole run, so what travels has to be summable.
///
/// # The sum is an integer, and it has to be
///
/// [`GenericAccumulators::merge`](super::accumulators::GenericAccumulators::merge) is
/// **order-independent** — its own test merges three shards in all six orders and asserts the
/// results agree — and that holds for its tables because they sum integer counts. A running `f64`
/// would have quietly broken it: floating-point addition is not associative, so a run that merged
/// its shards in a different order would produce a denominator differing in the last bits, and the
/// whole run's genotypes with it. Nothing would have failed, because no existing test compares this
/// field.
///
/// So the log error is counted in **fixed point, in units of 2⁻²⁰ ≈ 9.5 × 10⁻⁷**, where addition
/// is exact and order cannot matter. **The cost is bounded and it is bounded on the number that is
/// used**: each conversion rounds by at most half a unit, and every conversion that rounds at all
/// has at least one read behind it — an observation no read is behind carries a `q_sum` of exactly
/// zero, which converts exactly — so however many loci and shards a run has, the error on the
/// *mean* stays under 2⁻²¹ ≈ 4.8 × 10⁻⁷, against a mean log error of order 5 to 20.
///
/// **`i128` rather than `i64`, and the eight bytes are free.** One read group of one human sample
/// reaches about 2.2 × 10¹⁶ scaled units — a billion reads at up to 21.4 nats each — which an
/// `i64` holds four hundred times over. But [`add`](Self::add) folds across *samples* as well as
/// shards, because a read group is a library and a library can hold more than one plant; past
/// about four hundred human-scale samples in one read group an `i64` would saturate, and
/// `saturating_add` pins rather than panicking, so the symptom would be a mean that was merely
/// wrong. `i128` puts that beyond any run — 7 × 10²¹ such samples — and the map holds one entry
/// per read group, so the width costs bytes and not megabytes.
#[derive(Copy, Clone, Default, PartialEq, Eq, Debug)]
pub struct MintedReadErrors {
    /// Σ over the reads of `ln P(this read is wrong)`, in units of 2⁻²⁰.
    ///
    /// **Zero is a real value and not an absence.** A read at Phred 0 contributes exactly zero,
    /// because `ln 1 = 0`, and the mate-overlap rule silences a losing mate by giving it exactly
    /// that quality. So a group whose sum is zero is not a group nobody measured; `reads` is what
    /// tells those apart.
    log_error_sum_scaled: i128,
    /// How many reads that sum ran over.
    pub reads: u64,
}

impl MintedReadErrors {
    /// The totals for one observation: its own `q_sum` over its own `num_obs` reads.
    #[must_use]
    pub fn of_observation(q_sum: f64, num_obs: u32) -> Self {
        Self {
            // **Rust's float-to-integer `as` saturates rather than wrapping, and sends `NaN` to
            // zero.** So an out-of-range `q_sum` shows up as an absurd mean rather than a plausible
            // one, but a `NaN` shows up as nothing at all — the quieter of the two, and the one
            // worth knowing about. Neither is reachable from the fold, which sums table lookups;
            // it is said because the field it reads is public.
            log_error_sum_scaled: (q_sum * PARTS_OF_ONE).round() as i128,
            reads: u64::from(num_obs),
        }
    }

    /// The mean minted error **probability** over these reads — the scale's denominator.
    ///
    /// `None` where no read was seen, which is the honest answer: a scale needs a denominator and
    /// there is none. A caller that treated it as one would divide the fitted rate by nothing.
    #[must_use]
    pub fn mean_error_probability(self) -> Option<f64> {
        self.mean_log_error().map(f64::exp)
    }

    /// The mean minted error in log space, which is what the mean above is the exponential of.
    ///
    /// Offered beside it because the caller works in log space and the round trip through `exp`
    /// and back costs a few units in the last place — the same reason §3.3's identity is stated to
    /// a tolerance rather than bitwise.
    #[must_use]
    pub fn mean_log_error(self) -> Option<f64> {
        (self.reads > 0)
            .then(|| self.log_error_sum_scaled as f64 / PARTS_OF_ONE / self.reads as f64)
    }

    /// Fold another shard's or another sample's totals for the same read group into these.
    ///
    /// **A read group crosses both**, so the run's denominator is every one of them added up: two
    /// plants sequenced in one library were prepared by one chemistry, which is the whole reason
    /// the scale is per read group rather than per sample. Exact, and therefore independent of the
    /// order the folds happen in.
    pub fn add(&mut self, other: Self) {
        self.log_error_sum_scaled = self
            .log_error_sum_scaled
            .saturating_add(other.log_error_sum_scaled);
        self.reads = self.reads.saturating_add(other.reads);
    }
}

/// Total one locus's minted error per read group, into `out`.
///
/// `out` is scratch — cleared, then filled — for the reason
/// [`count_by_read_group`](super::depth_and_alt_reads::count_by_read_group)'s is: this runs once
/// per covered position over hundreds of millions of them.
///
/// **The same observations the read-group histogram counts, under the same gate.** A non-generic
/// locus contributes nothing here because it contributes nothing there, and the reads are the
/// complete witnesses for the same reason: a partial read's `q_sum` is over the stretch it saw,
/// and the fitted rate was never fitted from one.
pub fn minted_error_by_read_group(
    locus: &SampleLocusObservations,
    out: &mut Vec<(ReadGroupId, MintedReadErrors)>,
) {
    out.clear();
    if !matches!(locus.kind, LocusKind::Generic) {
        return;
    }
    for observation in locus.complete_observations() {
        let group = observation.read_group;
        let at = match out.iter().position(|&(seen, _)| seen == group) {
            Some(at) => at,
            None => {
                out.push((group, MintedReadErrors::default()));
                out.len() - 1
            }
        };
        out[at].1.add(MintedReadErrors::of_observation(
            observation.q_sum,
            observation.num_obs,
        ));
    }
    out.sort_unstable_by_key(|&(group, _)| group);
}

/// Fold one locus's totals into a running per-read-group table.
///
/// **Ascending read-group order, and the sum runs in it.** `f64` addition is not associative, so
/// two runs that visited a locus's groups in different orders would produce denominators that
/// differ in the last few bits, and the run's output would not be reproducible
/// (`doc/devel/ng/spec/read_likelihoods.md` §8). A `BTreeMap` gives that order for free and
/// [`minted_error_by_read_group`] sorts each locus's entries before they arrive.
pub fn fold_into(
    running: &mut BTreeMap<ReadGroupId, MintedReadErrors>,
    of_locus: &[(ReadGroupId, MintedReadErrors)],
) {
    for &(group, totals) in of_locus {
        running.entry(group).or_default().add(totals);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::locus_generation::{LocusLen, ReadWitness, SequenceObservation};
    use crate::ng::types::{ContigId, GenomeRegion, Position};

    fn observation(read_group: u32, num_obs: u32, q_sum: f64) -> SequenceObservation {
        SequenceObservation {
            bases: b"A"[..].into(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(read_group),
            num_obs,
            num_fwd: 0,
            q_sum,
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn locus(kind: LocusKind, observations: Vec<SequenceObservation>) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(10),
                end: Position(10),
            },
            reference_bases: b"A"[..].into(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind,
        }
    }

    /// **The mean is the mean log error, and the two numbers behind it are per read group.**
    ///
    /// Two libraries at one site, with different qualities and different depths, so neither the
    /// sum nor the count can be read off the other and pooling them would be visible.
    #[test]
    fn each_read_groups_mean_is_over_its_own_reads() {
        let mut out = Vec::new();
        minted_error_by_read_group(
            &locus(
                LocusKind::Generic,
                vec![
                    observation(0, 4, -12.0),
                    observation(1, 2, -3.0),
                    observation(0, 1, -6.0),
                ],
            ),
            &mut out,
        );

        assert_eq!(
            out.iter()
                .map(|&(group, _)| group.get())
                .collect::<Vec<_>>(),
            vec![0, 1],
            "ascending read-group order, which is what makes the fold reproducible",
        );
        // Group 0: two observations, 5 reads, −18 in total, so −3.6 a read.
        assert_eq!(out[0].1.reads, 5);
        assert_eq!(out[0].1.mean_log_error(), Some(-3.6));
        // Group 1: 2 reads at −3 in total, so −1.5 a read — and `exp(−1.5)` is what the
        // scale divides into the fitted rate.
        assert_eq!(out[1].1.reads, 2);
        assert_eq!(out[1].1.mean_log_error(), Some(-1.5));
        assert!(
            (out[1].1.mean_error_probability().expect("two reads") - (-1.5_f64).exp()).abs()
                < 1e-12,
        );
    }

    /// **A read group with no read has no mean**, which is the difference between a library that
    /// was measured and one that was not: a scale wanting a denominator gets no number rather
    /// than a one.
    #[test]
    fn no_reads_is_no_mean_rather_than_zero() {
        let empty = MintedReadErrors::default();
        assert_eq!(empty.mean_log_error(), None);
        assert_eq!(empty.mean_error_probability(), None);

        // And a sum of exactly zero **with** reads behind it is a real answer: a base the
        // instrument disclaims is charged an error probability of one, and `ln 1 = 0`.
        let disclaimed = MintedReadErrors::of_observation(0.0, 3);
        assert_eq!(disclaimed.mean_log_error(), Some(0.0));
        assert_eq!(disclaimed.mean_error_probability(), Some(1.0));
    }

    /// **The fold is exact and therefore order-independent**, which is the property
    /// [`GenericAccumulators::merge`](super::super::accumulators::GenericAccumulators) already
    /// promises for its tables and which an `f64` running sum here would have broken silently.
    ///
    /// The values are deliberately ones whose `f64` sum is order-dependent: `0.1`, `0.2` and
    /// `0.3` added left to right and right to left differ in the last bit. Here every order
    /// gives one number, bit for bit.
    #[test]
    fn folding_in_any_order_gives_the_same_totals() {
        let parts = [
            MintedReadErrors::of_observation(-0.1, 1),
            MintedReadErrors::of_observation(-0.2, 1),
            MintedReadErrors::of_observation(-0.3, 1),
            MintedReadErrors::of_observation(-1e-6, 1),
        ];
        // The premise: these four in two orders really do disagree as plain `f64` sums.
        let forwards: f64 = [-0.1, -0.2, -0.3, -1e-6].iter().sum();
        let backwards: f64 = [-1e-6, -0.3, -0.2, -0.1].iter().sum();
        assert_ne!(
            forwards, backwards,
            "if these agreed as floats the test would prove nothing",
        );

        let mut seen = Vec::new();
        for order in [
            [0, 1, 2, 3],
            [3, 2, 1, 0],
            [1, 3, 0, 2],
            [2, 0, 3, 1],
            [0, 3, 1, 2],
        ] {
            let mut total = MintedReadErrors::default();
            for at in order {
                total.add(parts[at]);
            }
            seen.push(total);
        }
        assert!(
            seen.windows(2).all(|pair| pair[0] == pair[1]),
            "the same four contributions in five orders: {seen:?}",
        );
        assert_eq!(seen[0].reads, 4);
    }

    /// **An observation no read is behind contributes nothing to either number**, which is what
    /// makes the rounding bound in [`MintedReadErrors`]'s doc unconditional: a conversion with no
    /// read behind it would add rounding error to a sum whose denominator did not grow.
    ///
    /// **It is not skipped, deliberately.** The error-rate histogram this is the denominator for
    /// does not skip one either, and the two paths visiting exactly the same observations is worth
    /// more than a branch that changes no number — the fold gives a read-less observation a `q_sum`
    /// of exactly zero, so it contributes zero and rounds exactly.
    #[test]
    fn an_observation_no_read_is_behind_contributes_nothing() {
        let mut out = Vec::new();
        minted_error_by_read_group(
            &locus(
                LocusKind::Generic,
                vec![observation(0, 0, 0.0), observation(0, 3, -6.0)],
            ),
            &mut out,
        );

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1.reads, 3, "the three reads that exist");
        assert_eq!(
            out[0].1.mean_log_error(),
            Some(-2.0),
            "and the mean is over those three and nothing else",
        );
    }

    /// **A repeat-tract locus contributes nothing**, because it contributes nothing to the table
    /// whose rate this is the denominator for — the STR path has its own noise model.
    #[test]
    fn only_generic_loci_are_counted() {
        let mut out = Vec::new();
        for kind in [LocusKind::SsrBundle, LocusKind::Generic] {
            let generic = matches!(kind, LocusKind::Generic);
            let named = format!("{kind:?}");
            minted_error_by_read_group(&locus(kind, vec![observation(0, 4, -12.0)]), &mut out);
            let counted = out.len();
            assert_eq!(
                counted,
                usize::from(generic),
                "{named} contributed {counted} read groups",
            );
        }
    }

    /// **A partial read is not counted**, for the same reason: its `q_sum` is over the stretch it
    /// saw, and the rate this is the denominator for was never fitted from one.
    #[test]
    fn a_read_that_ran_out_is_not_counted() {
        let mut ran_out = observation(0, 4, -12.0);
        ran_out.read_witness = ReadWitness::from_left(
            1,
            LocusLen::of_region(GenomeRegion {
                contig: ContigId(0),
                start: Position(10),
                end: Position(10),
            }),
        )
        .expect("a one-base run inside a one-base locus");

        let mut out = Vec::new();
        minted_error_by_read_group(
            &locus(LocusKind::Generic, vec![ran_out, observation(0, 2, -4.0)]),
            &mut out,
        );

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1.reads, 2, "only the complete observation's reads");
        assert_eq!(out[0].1.mean_log_error(), Some(-2.0));
    }

    /// **The fixed-point sum's error on the mean stays under half a unit however much is added**,
    /// which is the bound the type's doc claims. Measured rather than asserted: a hundred
    /// thousand reads of an awkward quality, against the exact answer.
    #[test]
    fn the_fixed_point_rounding_does_not_accumulate_in_the_mean() {
        let awkward = -7.123_456_789_012_345;
        let mut total = MintedReadErrors::default();
        for _ in 0..100_000 {
            total.add(MintedReadErrors::of_observation(awkward, 1));
        }
        let error = (total.mean_log_error().expect("reads") - awkward).abs();
        assert!(
            error < 1.0 / PARTS_OF_ONE / 2.0,
            "the mean is off by {error}, past half a unit of 2^-20",
        );
        assert_eq!(total.reads, 100_000);
    }
}
