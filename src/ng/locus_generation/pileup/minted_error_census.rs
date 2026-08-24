//! **Measurement scaffolding, not part of the walk** — per-read-group sums of the minted
//! per-read error, kept in *both* the shapes an average can be taken in.
//!
//! The pre-pass's calibration accumulator
//! ([`calibration`](crate::ng::parameter_estimation::generic::calibration)) keeps Σ `ln ε` and the
//! read count, and the scale it feeds divides by their **geometric** mean, `exp(Σ ln ε / n)`. The
//! **arithmetic** mean, `Σ ε / n`, is not recoverable from those two numbers and nothing in the
//! walk carries it, so the question *how far apart are the two on real reads* could only ever be
//! answered by summing `ε` where a read still exists — which is here, and nowhere later.
//!
//! # What it counts, and why it is the same reads the accumulator sees
//!
//! `OpenPileupRecord::finalise` resolves every folded read's witness against
//! the record's final footprint on its way to building the observations. At that point a read's
//! own `contribution.q_sum` is exactly its own `ln ε` — one read, one contribution, carried across
//! widens untouched — and its witness and read group are settled. Recording there, for
//! [`ReadWitness::Complete`](crate::ng::locus_generation::ReadWitness) reads only, gives exactly
//! the reads that
//! [`minted_error_by_read_group`](crate::ng::parameter_estimation::generic::calibration::minted_error_by_read_group)
//! later sums over: complete witnesses at generic loci, before the per-position depth cap, which
//! draws on counts and never on a quality.
//!
//! **That identity is checked rather than asserted.** `examples/ng_minted_error_means.rs` prints
//! this census's Σ `ln ε` and read count beside the accumulator's own, and a run where they
//! disagree says so in its output.
//!
//! # A locus is counted when it is kept, not when it is built
//!
//! Two things would otherwise break that identity, and both were found by the check above rather
//! than reasoned about. **The walk has two paths that emit observations** — the general fold
//! through `open_record.rs`'s `OpenPileupRecord` and the ordinary-column fast lane in
//! `fast_column.rs`, which never builds a record at all — so both call
//! [`record_read`]. And **the walk builds records in a halo beyond the region it was asked
//! for**, which the generator then discards by their start position
//! (`PileupGeneratorCounts::records_outside_region`). So a read's totals land in a **pending**
//! bucket keyed by its locus's start and move into the table only when the generator keeps that
//! locus. Committing at build time instead overcounted HG002's 100 benchmark regions at 300× by
//! 1,489,219 reads in 174 million.
//!
//! # Off by default, and it costs one atomic load per record when it is
//!
//! Armed with `PVC_MINTED_ERROR_CENSUS=1`, read once per process. A record's finalise checks
//! [`enabled`] once, not once per read, so an unarmed walk pays a cached load per record against
//! the ~90 reads that record folds. Armed, it takes a mutex per read, which is why this is a
//! measurement mode and not a counter.

use std::collections::BTreeMap;
use std::sync::{Mutex, OnceLock};

use crate::ng::types::ReadGroupId;

/// One read group's minted error in both summable shapes, over the reads this census saw.
///
/// **Both sums are `f64` and their low bits are order-dependent**, which is exactly what the
/// pre-pass's own accumulator refuses — it holds fixed point so that merging region shards in any
/// order gives one answer. Here that does not matter and paying for it would: this census runs in
/// one thread over one sample and its answer is quoted to three significant figures.
#[derive(Copy, Clone, Default, PartialEq, Debug)]
pub struct MintedErrorTotals {
    /// Σ over the reads of `ln P(this read is wrong)` — the same quantity an observation's
    /// `q_sum` sums, before it is pooled by allele.
    pub log_error_sum: f64,
    /// Σ over the reads of `P(this read is wrong)` itself. **This is the number that exists
    /// nowhere else**: the fold discards the individual reads, so `Σ ε` cannot be recovered from
    /// `Σ ln ε` afterwards.
    pub error_sum: f64,
    /// How many reads both sums ran over.
    pub reads: u64,
    /// Of those, how many carry `ε` of **exactly one** — a read charged a full unit of error.
    ///
    /// **Kept apart because these dominate the arithmetic mean and barely touch the geometric
    /// one**, so the gap between the two means cannot be read without them. A read reaches
    /// `ε = 1` when its window's base quality is Phred 0, which is what the mate-overlap rule
    /// gives the losing mate of an overlapping pair: it is silenced, not dropped, so it still
    /// counts as a read. `ln 1 = 0`, so such a read adds nothing to the log sum and one whole
    /// unit to the probability sum.
    pub reads_charged_a_full_unit: u64,
}

impl MintedErrorTotals {
    /// The geometric mean of the per-read error probabilities — what the scale's denominator is
    /// (`doc/devel/ng/spec/read_likelihoods.md` §3.2). `None` where no read was seen.
    #[must_use]
    pub fn geometric_mean(self) -> Option<f64> {
        (self.reads > 0).then(|| (self.log_error_sum / self.reads as f64).exp())
    }

    /// The arithmetic mean of the same probabilities — what §3.2 asked for before the correction
    /// of 2026-08-24. `None` where no read was seen.
    #[must_use]
    pub fn arithmetic_mean(self) -> Option<f64> {
        (self.reads > 0).then(|| self.error_sum / self.reads as f64)
    }

    /// **How much of the arithmetic mean is reads charged a full unit of error** — their share of
    /// `Σ ε`, which is their count over that sum because each contributes exactly one.
    ///
    /// The number that says whether the arithmetic mean is measuring the chemistry or measuring
    /// how often a mate pair overlaps. `None` where no read was seen or nothing was charged.
    #[must_use]
    pub fn full_unit_share_of_arithmetic_mean(self) -> Option<f64> {
        (self.reads > 0 && self.error_sum > 0.0)
            .then(|| self.reads_charged_a_full_unit as f64 / self.error_sum)
    }
}

/// Whether the census is armed for this process. Read once; an unarmed walk pays a cached load
/// per record.
#[must_use]
pub fn enabled() -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| std::env::var_os("PVC_MINTED_ERROR_CENSUS").is_some_and(|v| v == "1"))
}

/// What the census holds: the loci the generator has kept, and the ones still waiting to hear.
///
/// A `Mutex` over both and not an atomic pair, because two `f64` sums have to move together for
/// a read: an interleaving that added `ln ε` from one read and `ε` from another would give the
/// two means different denominators without changing the count.
#[derive(Default)]
struct Census {
    /// The answer, per read group, over the loci the generator kept.
    kept: BTreeMap<ReadGroupId, MintedErrorTotals>,
    /// Built but not yet ruled on, keyed by the locus's start position. The walker buffers
    /// finished loci in a queue, so several can be pending at once and a plain "the current
    /// locus" buffer would attribute one locus's reads to another.
    pending: BTreeMap<u64, BTreeMap<ReadGroupId, MintedErrorTotals>>,
}

fn census() -> &'static Mutex<Census> {
    static CENSUS: OnceLock<Mutex<Census>> = OnceLock::new();
    CENSUS.get_or_init(|| Mutex::new(Census::default()))
}

/// Record one read's minted error at the locus starting at `locus_start`.
///
/// `log_error` is the read's own `ln ε` — one read, one contribution, before anything pools it.
/// The caller has already checked [`enabled`] and that the read's witness is complete; this does
/// not re-check either, so that an unarmed walk never reaches here at all.
pub fn record_read(locus_start: u64, read_group: ReadGroupId, log_error: f64) {
    let mut census = census()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let totals = census
        .pending
        .entry(locus_start)
        .or_default()
        .entry(read_group)
        .or_default();
    totals.log_error_sum += log_error;
    totals.error_sum += log_error.exp();
    totals.reads += 1;
    // `>= 0.0` and not `== 0.0`: `phred_to_ln_perr(0)` is pinned to `+0.0` and the mint takes a
    // `max`, so a silenced read arrives as exactly `+0.0` and nothing can arrive above it — the
    // comparison is written this way so that a positive value, which would mean a probability
    // above one, lands in this count rather than escaping it.
    if log_error >= 0.0 {
        totals.reads_charged_a_full_unit += 1;
    }
}

/// The generator kept the locus starting at `locus_start`: move its reads into the answer.
///
/// A start with nothing pending is a locus no read reached, and is not an error.
pub fn keep_locus(locus_start: u64) {
    let mut census = census()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let Some(of_locus) = census.pending.remove(&locus_start) else {
        return;
    };
    for (group, totals) in of_locus {
        let running = census.kept.entry(group).or_default();
        running.log_error_sum += totals.log_error_sum;
        running.error_sum += totals.error_sum;
        running.reads += totals.reads;
        running.reads_charged_a_full_unit += totals.reads_charged_a_full_unit;
    }
}

/// The generator discarded the locus starting at `locus_start` — it was built in the halo beyond
/// the region asked for. Its reads leave without being counted.
pub fn drop_locus(locus_start: u64) {
    census()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .pending
        .remove(&locus_start);
}

/// How many loci were built and never ruled on — a walk that ended mid-region, or a locus the
/// generator neither kept nor discarded.
///
/// **Printed by the tool rather than asserted**, because a non-zero here does not make the
/// answer wrong: those reads are simply not in it. It is the number that says how much.
#[must_use]
pub fn loci_never_ruled_on() -> usize {
    census()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .pending
        .len()
}

/// The kept loci's totals, in ascending read-group order.
#[must_use]
pub fn snapshot() -> Vec<(ReadGroupId, MintedErrorTotals)> {
    census()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
        .kept
        .iter()
        .map(|(&group, &totals)| (group, totals))
        .collect()
}

/// Empty both tables, so one process can walk several samples and report each on its own.
pub fn reset() {
    let mut census = census()
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    census.kept.clear();
    census.pending.clear();
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The two means are different numbers, and the fixture is chosen so that a swap is
    /// visible.** One read at Phred 20 (`ε = 0.01`) and one at Phred 40 (`ε = 0.0001`) is
    /// `read_likelihoods.md` §3.2's own example: arithmetic 0.00505, geometric 0.001, a factor
    /// of 5.05 between them.
    #[test]
    fn the_two_means_differ_by_the_specs_own_example() {
        let mut totals = MintedErrorTotals::default();
        for error in [0.01_f64, 0.000_1] {
            totals.log_error_sum += error.ln();
            totals.error_sum += error;
            totals.reads += 1;
        }

        let arithmetic = totals.arithmetic_mean().expect("two reads");
        let geometric = totals.geometric_mean().expect("two reads");
        assert!(
            (arithmetic - 0.005_05).abs() < 1e-12,
            "arithmetic mean came back {arithmetic}",
        );
        assert!(
            (geometric - 0.001).abs() < 1e-12,
            "geometric mean came back {geometric}",
        );
        assert!(
            (arithmetic / geometric - 5.05).abs() < 1e-9,
            "the ratio came back {}",
            arithmetic / geometric,
        );
    }

    /// **No read is no mean, not a zero** — the same rule the pre-pass's accumulator follows,
    /// because a ratio of the two means needs both to exist.
    #[test]
    fn no_reads_is_no_mean() {
        let empty = MintedErrorTotals::default();
        assert_eq!(empty.geometric_mean(), None);
        assert_eq!(empty.arithmetic_mean(), None);
    }

    /// **One silenced read among many good ones is most of the arithmetic mean and almost none
    /// of the geometric one**, which is the whole reason the full-unit count is kept apart.
    ///
    /// Ninety-nine reads at `ε = 10⁻⁴` and one at `ε = 1`: the arithmetic mean is 1.0099 in a
    /// hundred and the silenced read is 99.0% of it, while the geometric mean is 1.096 in ten
    /// thousand — the silenced read moved it by a tenth, not by a hundredfold.
    #[test]
    fn one_silenced_read_owns_the_arithmetic_mean_and_barely_moves_the_geometric_one() {
        let mut totals = MintedErrorTotals::default();
        for _ in 0..99 {
            totals.log_error_sum += 0.000_1_f64.ln();
            totals.error_sum += 0.000_1;
            totals.reads += 1;
        }
        totals.error_sum += 1.0;
        totals.reads += 1;
        totals.reads_charged_a_full_unit += 1;

        let arithmetic = totals.arithmetic_mean().expect("a hundred reads");
        let geometric = totals.geometric_mean().expect("a hundred reads");
        assert!(
            (arithmetic - 0.010_099).abs() < 1e-9,
            "arithmetic mean came back {arithmetic}",
        );
        assert!(
            (geometric - 1.096_478_196e-4).abs() < 1e-12,
            "geometric mean came back {geometric}",
        );
        let share = totals
            .full_unit_share_of_arithmetic_mean()
            .expect("something was charged");
        assert!(
            (share - 0.990_197).abs() < 1e-6,
            "the silenced read's share came back {share}",
        );
        // And the ratio the tool prints: 92.1, from one read in a hundred.
        assert!(
            (arithmetic / geometric - 92.104).abs() < 1e-3,
            "the ratio came back {}",
            arithmetic / geometric,
        );
    }

    /// **Reads that all carry the same error make the two means agree exactly**, which is the
    /// null this measurement is against: the gap between them is a fact about the *spread* of
    /// the per-read errors and about nothing else.
    #[test]
    fn one_repeated_error_makes_the_two_means_equal() {
        let mut totals = MintedErrorTotals::default();
        for _ in 0..97 {
            totals.log_error_sum += 0.003_f64.ln();
            totals.error_sum += 0.003;
            totals.reads += 1;
        }
        let arithmetic = totals.arithmetic_mean().expect("reads");
        let geometric = totals.geometric_mean().expect("reads");
        assert!(
            (arithmetic / geometric - 1.0).abs() < 1e-12,
            "identical reads gave a ratio of {}",
            arithmetic / geometric,
        );
    }
}
