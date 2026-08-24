//! **The denominator of the read likelihood's error-rate scale, per read group.**
//!
//! **The caller holds one error probability per observation, not one per read.** The merge keeps,
//! for each allele in each read group at each locus, how many reads support it and the sum of their
//! log error probabilities; the reads are gone from there on. So the model charges
//! `exp(q_sum / num_obs)` — the geometric mean of those reads' minted errors — and the scale is one
//! addition in log space per observation, never a multiplication read by read
//! (`doc/devel/ng/spec/read_likelihoods.md` §3.2, §3.3). *(§3.2 states the rule per read, and may:
//! scaling every read and scaling their geometric mean are the same operation,
//! `exp(Σ ln(s·ε) / n) = s · exp(Σ ln ε / n)`.)*
//!
//! The scale makes that charged average come out at the rate the pre-pass measured for the library:
//!
//! ```text
//!                     fitted error rate for this read group
//! scale  =  ──────────────────────────────────────────────────────────
//!           geometric mean of the minted error over that group's reads
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
//! **It is not that the arithmetic mean would have been a worse choice — there was no place to use
//! it.** Nowhere in the model does a per-read `ε` survive to be averaged arithmetically; what the
//! model charges is `exp(q_sum / num_obs)`, and so does production
//! (`var_calling/posterior_engine.rs`, which has no recalibration at all — there is nothing to copy
//! there but the quantity). A scale built from an arithmetic mean and applied to a geometric one
//! would not make the calibrated property hold in the model's own terms, so supplying the
//! arithmetic sum would have bought an inexactness rather than removed one — at the price of a
//! second accumulation at fold time and a field on every observation. Corrected in the spec on
//! 2026-08-24 by the owner.
//!
//! # The sites are the same; the reads are not, and the gap is 3 parts in 100 at 300×
//!
//! §3.2 requires the average to run over exactly the reads the fitted rate was fitted from. **The
//! *sites* are identical by construction** — both paths run behind one `LocusKind::Generic` gate in
//! [`add_locus`](super::accumulators::GenericAccumulators::add_locus), before its inbreeding-mode
//! branch, and both iterate `complete_observations()`. Neither the library count, nor the ploidy
//! map, nor a supplied inbreeding coefficient divides them.
//!
//! **The per-position depth cap does.** The histogram route thins every position to at most
//! [`MAX_BINNED_DEPTH`](super::depth_bins) reads before fitting; this fold thins nothing. Per site
//! that is harmless — the draw is hypergeometric on counts and never looks at a read's quality —
//! but across sites it re-weights: a 500-read position casts 500 votes here and 124 in the
//! population the rate was fitted from, and deep positions are where reads pile up from
//! elsewhere and where mapping quality collapses.
//!
//! **Measured, not argued** (`examples/ng_minted_error_means.rs`). On HG002's 100 benchmark regions
//! at 300×, where the fit sees 41 read-positions in 100, the denominator's geometric mean is
//! 2.9055 × 10⁻⁴ against 2.9862 × 10⁻⁴ with each position thinned first — **2.7%, or 0.12 Phred**.
//! On the tomato cohort — 2.5× to 28.6× over 63 accessions — it is nothing: on the deepest of them
//! 228,468,065 read-positions of 228,492,796 are under the cap and the mean moves by 1.0000.
//! **This fold does not thin, and that is decided rather than pending** (owner, 2026-08-24): the
//! scale is applied to every read at calling time, so the average it is built from is over every
//! read. The 2.7% is carried knowingly, and it is a question about how the *fit* weights deep sites
//! against shallow ones rather than about this average (spec §3.2).
//!
//! # One thing the numerator can be that this cannot
//!
//! A read group standing on fewer than [`MIN_SITES_TO_FIT`](super::MIN_SITES_TO_FIT) sites — ten
//! thousand — does not get its own fitted rate:
//! [`resolve_error_rates`](super::fallback::resolve_error_rates) hands it the mean of the other
//! groups' rates, or a supplied one, or a default. **Its denominator is still its own reads**, so for such a group the
//! scale is "make this library's average charged error come out at somebody else's measured rate".
//! That may be exactly right — it is what borrowing a rate means — but §3.2's sentence about one
//! site set does not describe it, and a capture panel or a minor library in a multi-library sample
//! reaches it.

use std::collections::BTreeMap;

use crate::ng::locus_generation::{LocusKind, SampleLocusObservations};
use crate::ng::types::ReadGroupId;

/// How many parts of one the log-error sum is counted in — see [`MintedReadErrors`].
const PARTS_OF_ONE: f64 = (1_u64 << 20) as f64;

/// One read group's minted error, summed over the sites a fit read from.
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
/// is exact and order cannot matter. **The cost is bounded, and it is bounded on the number that
/// is used**: each conversion rounds by at most half a unit, and every conversion the *walk*
/// produces has at least one read behind it — the fold creates an observation only when a read
/// arrives and then counts that read — so however many loci and shards a run has, the error on
/// the *mean* is **at most** 2⁻²¹ ≈ 4.768 × 10⁻⁷, against a mean log error of order 5 to 20.
///
/// **At most, and the bound is attained rather than approached.** A `q_sum` of exactly half a
/// unit — `−(2k+1)·2⁻²¹` — rounds the whole half unit away from zero at every conversion, and the
/// miss on the mean is then 4.768 × 10⁻⁷ at one observation and still 4.768 × 10⁻⁷ at twenty
/// million: it does not grow with the run, and it does not shrink. Reads sharing a conversion
/// make it better, never worse, because there are fewer conversions per read.
///
/// **The bound is about the walk's observations, not about every value this type will accept.**
/// [`of_observation`](Self::of_observation) is public and will take a `q_sum` with `num_obs = 0`
/// — rounding a sum that no read enlarges the denominator for — and a hand-written fixture can
/// hand it one. That is not reachable from the fold on either column path, and stating it here
/// rather than defending against it is deliberate: a guard would turn a fixture's mistake into
/// silence.
///
/// **`i128` rather than `i64`, and an `i64` would not have been enough for one deep human
/// sample.** What [`reads`](Self::reads) counts is a read **at a position**, not a read — an
/// observation contributes its `num_obs` at every locus it appears at — so a sample's total is
/// its covered length times its depth, not its read count. Measured on HG002's 100 benchmark
/// regions at 300×: 172,616,054 over 571,984 bases, which is 301.8 a base
/// (`examples/ng_minted_error_means.rs`). So:
///
/// **The largest one of these can get is one sample's whole genome**, because a read group belongs
/// to exactly one sample — `build_read_groups` mints one entry per declared read group per file and
/// `group_by_sample` files each under the single sample its header names, so an identifier cannot
/// span samples. So the bound to clear is a single deep library, not a cohort:
///
/// - a human genome at 30× is about 9.3 × 10¹⁰ read-positions, and at the mean log error measured
///   on that sample (8.145 nats) that is 7.9 × 10¹⁷ scaled units — an `i64` has 11.6× headroom;
/// - **the same genome at 300× is 7.9 × 10¹⁸, which is 86% of `i64::MAX` — 1.16× headroom, for one
///   library.** A deeper run, a larger genome or a noisier library goes over.
///
/// **A single read's contribution is bounded by the *smaller* of its two qualities, not the
/// larger**, which is easy to get backwards: the mint takes `max(ln ε_BQ, ln ε_MQ)` of two
/// negative numbers, so it returns the one nearer zero. Over every `(base quality, mapping
/// quality)` byte pair the most negative value the mint can return is −58.716 — `phred_to_ln_perr`
/// at 255 — and on BWA output, where mapping quality tops out at 60, the real ceiling is
/// **−13.816**. **At that ceiling one 300× human library needs 1.35 × 10¹⁹ scaled units, which an
/// `i64` does not hold at all.**
///
/// So `i64` is not comfortably enough for a single deep library, and `saturating_add` pins rather
/// than panicking — the symptom would have been a mean that was merely wrong rather than a run that
/// stopped. `i128` puts it beyond any depth. It costs 16 bytes a read group and not 8, because
/// `i128` also aligns the struct to 16 — 32 bytes against `i64`'s 16 — and the map holds one entry
/// per read group.
///
/// [`merge`]: super::accumulators::GenericAccumulators::merge
#[derive(Copy, Clone, Default, PartialEq, Eq, Debug)]
pub struct MintedReadErrors {
    /// Σ over the reads of `ln P(this read is wrong)`, in units of 2⁻²⁰.
    ///
    /// **Zero is a real value and not an absence.** A read at Phred 0 contributes exactly zero,
    /// because `ln 1 = 0`, and the mate-overlap rule silences a losing mate by giving it exactly
    /// that quality. So a group whose sum is zero is not a group nobody measured; `reads` is what
    /// tells those apart.
    log_error_sum_scaled: i128,
    /// How many reads that sum ran over — **a read at a position, counted once for every
    /// position it is seen at**, which is what an observation's `num_obs` is and what the fitted
    /// error rate this is the denominator for is a rate *per*. A 150-base read covering 150
    /// generic loci is 150 of these. So a sample's total is its covered length times its depth:
    /// 172,616,054 over 571,984 bases on HG002 at 300×.
    reads: u64,
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

    /// How many reads this ran over — see the field for what a "read" is counted as here.
    ///
    /// **Read-only, like the sum beside it**, and that is the point: the mean is meaningful only
    /// because the two numbers cover the same reads, and a public field would let a holder move
    /// one without the other. Both change together or not at all, through
    /// [`of_observation`](Self::of_observation) and [`add`](Self::add).
    #[must_use]
    pub fn reads(self) -> u64 {
        self.reads
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
    /// and back costs a few units in the last place — the same reason §3.6's identity with §3.3
    /// is stated to a tolerance rather than bitwise. (§3.3's own aggregation identity *is*
    /// bitwise; it is the contamination mixture's that is not.)
    #[must_use]
    pub fn mean_log_error(self) -> Option<f64> {
        (self.reads > 0)
            .then(|| self.log_error_sum_scaled as f64 / PARTS_OF_ONE / self.reads as f64)
    }

    /// Fold another region shard's totals for the same read group into these.
    ///
    /// Exact, and therefore independent of the order the folds happen in.
    ///
    /// **Shards of one sample, and there is nothing else to fold.** A read group belongs to
    /// exactly one sample: `build_read_groups` mints one entry per declared read group per file
    /// and `group_by_sample` files each identifier under the single sample its header names, so
    /// two samples cannot share one. The scale is per read group rather than per sample because
    /// chemistry belongs to the library — a sample sequenced twice has two of these — not because
    /// a library spans samples.
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
/// **The order the groups arrive in cannot change the answer, and that is the fixed point's
/// doing rather than this function's.** Each group is its own key and each sum is an exact
/// integer, so nothing here is order-sensitive — measured, not argued: deleting
/// [`minted_error_by_read_group`]'s sort changes no total in either this module's tests or
/// [`accumulators`](super::accumulators)'. The sort is kept because that function *states* it
/// returns ascending order and a caller reading the scratch vector directly should get it —
/// `each_read_groups_mean_is_over_its_own_reads` hands it a locus whose groups arrive descending,
/// so the claim is pinned. **It is not what makes this fold reproducible.** Were the sum an
/// `f64` it would have to be, because `f64` addition is not associative — which is the trap
/// [`MintedReadErrors`] exists to sidestep (`doc/devel/ng/spec/read_likelihoods.md` §8 on
/// determinism).
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
    use crate::ng::locus_generation::{LocusLen, ReadWitness, SequenceObservation, SsrDetail};
    use crate::ng::types::{ContigId, GenomeRegion, Motif, Position};

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
    ///
    /// **The higher-numbered group is written first on purpose.** With the groups arriving
    /// ascending, the sort at the end of [`minted_error_by_read_group`] is a no-op and the order
    /// assertion below passes whether the sort is there or not — which is what it was, until a
    /// mutation showed that deleting the sort left every test green.
    #[test]
    fn each_read_groups_mean_is_over_its_own_reads() {
        let mut out = Vec::new();
        minted_error_by_read_group(
            &locus(
                LocusKind::Generic,
                vec![
                    observation(1, 2, -3.0),
                    observation(0, 4, -12.0),
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
            "ascending read-group order, which this function states it returns",
        );
        // Group 0: two observations, 5 reads, −18 in total, so −3.6 a read.
        assert_eq!(out[0].1.reads(), 5);
        assert_eq!(out[0].1.mean_log_error(), Some(-3.6));
        // Group 1: 2 reads at −3 in total, so −1.5 a read — and `exp(−1.5)` is what the
        // scale divides into the fitted rate.
        assert_eq!(out[1].1.reads(), 2);
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
        assert_eq!(seen[0].reads(), 4);
    }

    /// **An observation no read is behind is not skipped, and what it costs depends on its
    /// `q_sum`.** The error-rate histogram this is the denominator for does not skip one either,
    /// and the two paths visiting exactly the same observations is worth more than a branch.
    ///
    /// **The walk cannot build one at all**, which is the fact the rounding bound in
    /// [`MintedReadErrors`]'s doc rests on: both column paths create an observation only when a
    /// read arrives and count that read in the same step. So the first case below — `q_sum` of
    /// zero, which converts exactly — is the only shape the fold produces. **The second case is
    /// the one a fixture can produce and the walk cannot**, and it is pinned rather than
    /// defended against, so that the cost is written down: a read-less observation carrying a
    /// non-zero `q_sum` moves the sum without moving the count, and the mean moves with it.
    #[test]
    fn an_observation_no_read_is_behind_is_counted_as_its_q_sum_says() {
        let mut out = Vec::new();
        minted_error_by_read_group(
            &locus(
                LocusKind::Generic,
                vec![observation(0, 0, 0.0), observation(0, 3, -6.0)],
            ),
            &mut out,
        );

        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1.reads(), 3, "the three reads that exist");
        assert_eq!(
            out[0].1.mean_log_error(),
            Some(-2.0),
            "and the mean is over those three and nothing else",
        );

        // The shape only a hand-written fixture reaches: no read, and a `q_sum` anyway. −1.5
        // lands on the sum, 0 lands on the count, and the mean over the three real reads moves
        // from −2.0 to −2.5. **That is the answer, not a bug to be caught here** — a guard would
        // turn a fixture's mistake into silence, and the walk does not produce this. (−1.5 and
        // −6.0 are both whole multiples of 2⁻²⁰, so the equality is exact rather than
        // approximate and a wrong sum cannot hide inside a tolerance.)
        minted_error_by_read_group(
            &locus(
                LocusKind::Generic,
                vec![observation(0, 0, -1.5), observation(0, 3, -6.0)],
            ),
            &mut out,
        );
        assert_eq!(
            out[0].1.reads(),
            3,
            "still three reads, because none arrived"
        );
        assert_eq!(
            out[0].1.mean_log_error(),
            Some(-2.5),
            "and the read-less observation's own q_sum is in the sum regardless",
        );
    }

    /// **A repeat-tract locus contributes nothing**, because it contributes nothing to the table
    /// whose rate this is the denominator for — the STR path has its own noise model.
    #[test]
    fn only_generic_loci_are_counted() {
        let mut out = Vec::new();
        // **The generic locus goes first**, so the scratch buffer is non-empty when the
        // repeat-tract locus arrives. With the empty one first, a build that emptied `out` only
        // *after* the kind gate — leaving the previous locus's totals standing — passed this and
        // every other test in both modules.
        for kind in [LocusKind::Generic, LocusKind::SsrBundle] {
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

    /// **The sum is wide enough for a run, and an `i64` would not have been.** Swapping
    /// `log_error_sum_scaled` to `i64` left every other test in this module and in
    /// [`accumulators`](super::accumulators) green, because their fixtures peak near 2.7 × 10⁷
    /// scaled units and an `i64` holds 9.2 × 10¹⁸. This is the test that sees it.
    ///
    /// The scale is a real one. `reads` counts a read **at a position**, so a human genome at 30×
    /// is about 9.3 × 10¹⁰ of them; the 1.288 × 10¹² accumulated here is about fourteen such
    /// samples in one read group, or one sample at 415×. At an average of 8 nats a read-position
    /// that is 1.08 × 10¹⁹ scaled units, past `i64::MAX`.
    ///
    /// **Every value is a whole power of two**, so the mean is exactly −8 and the assertion needs
    /// no tolerance to hide behind: an `i64` `saturating_add` pins at `i64::MIN` and returns
    /// −6.826666…, which is not a number a tolerance could excuse.
    #[test]
    fn the_sum_is_wide_enough_for_a_run_that_an_i64_would_have_saturated() {
        // 2^31 read-positions at exactly 8 nats each: `q_sum` is −2^34, which scales to −2^54.
        let reads_per_call: u32 = 1 << 31;
        let q_sum = -8.0 * f64::from(reads_per_call);
        let mut total = MintedReadErrors::default();
        for _ in 0..600 {
            total.add(MintedReadErrors::of_observation(q_sum, reads_per_call));
        }

        assert_eq!(total.reads(), 600 * u64::from(reads_per_call));
        assert_eq!(
            total.mean_log_error(),
            Some(-8.0),
            "600 × 2^54 scaled units is 1.08e19, and i64::MAX is 9.22e18 — an i64 here pins at \
             i64::MIN and gives -6.826666666666666",
        );
    }

    /// **A repeat tract contributes nothing either**, which is the variant the test above does
    /// not name: [`LocusKind`] has three arms — `Generic`, `Ssr`, `SsrBundle` — and a gate written
    /// as "everything but a bundle" passes that test while pouring every repeat-tract read into
    /// the denominator of a rate fitted from no repeat-tract read at all.
    #[test]
    fn a_repeat_tract_contributes_nothing() {
        let mut out = Vec::new();
        let tract = LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a two-base motif"),
            left_flank: Box::from(&b"CCCGGG"[..]),
            right_flank: Box::from(&b"TTTAAA"[..]),
        });
        minted_error_by_read_group(&locus(tract, vec![observation(0, 4, -12.0)]), &mut out);
        assert_eq!(
            out.len(),
            0,
            "a repeat tract reached the error-rate denominator: {out:?}",
        );
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
        assert_eq!(out[0].1.reads(), 2, "only the complete observation's reads");
        assert_eq!(out[0].1.mean_log_error(), Some(-2.0));
    }

    /// **The fixed-point sum's error on the mean stays within half a unit however much is
    /// added**, which is the bound the type's doc claims. Measured rather than asserted: a
    /// hundred thousand reads of an awkward quality, against the exact answer.
    ///
    /// **The tolerance is written as a literal `2⁻²¹` on purpose.** Spelled `1.0 / PARTS_OF_ONE
    /// / 2.0` it moves with the constant it is checking, so coarsening the grid from 2⁻²⁰ to
    /// 2⁻¹⁹ — a mutation that really does change the answer, a `q_sum` of −17.223 620 414 733 887
    /// becoming −17.223 619 461 059 570 — left this test green. Written as a number, the grain
    /// is what the test pins.
    #[test]
    fn the_fixed_point_rounding_does_not_accumulate_in_the_mean() {
        let awkward = -7.123_456_789_012_345;
        let mut total = MintedReadErrors::default();
        for _ in 0..100_000 {
            total.add(MintedReadErrors::of_observation(awkward, 1));
        }
        let error = (total.mean_log_error().expect("reads") - awkward).abs();
        assert!(
            error <= 4.768_371_582_031_25e-7,
            "the mean is off by {error}, past half a unit of 2^-20",
        );
        assert_eq!(total.reads(), 100_000);

        // **The bound is reached, not merely respected**, which is why it is `<=` above. A
        // `q_sum` sitting exactly half a unit off the grid rounds the whole half unit away from
        // zero at every conversion, and one read is enough to see it.
        let exactly_half_a_unit = -3.0 / (2.0 * PARTS_OF_ONE);
        let worst = MintedReadErrors::of_observation(exactly_half_a_unit, 1);
        let miss = (worst.mean_log_error().expect("one read") - exactly_half_a_unit).abs();
        assert!(
            (miss - 4.768_371_582_031_25e-7).abs() < 1e-18,
            "the worst case missed by {miss}, which is not the bound",
        );
    }
}
