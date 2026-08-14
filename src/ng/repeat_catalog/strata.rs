//! Counting and sampling STR loci by stratum — what the parameter pre-pass needs and the
//! reason the catalog exists at all.
//!
//! A **stratum** is one `(period, repeat count)` pair. Fitting the stutter model needs `cap`
//! loci drawn at random from each of them, and the rule that gives exactly `cap` and a
//! uniform sample is "the `cap` lowest hashes in the stratum" — which needs the stratum
//! enumerated first (`parameter_prepass_joint_loci.md` §3.2). That enumeration is this file.

use std::collections::HashMap;

use crate::ng::region_typing::segment_criteria::SsrSegment;

/// How many loci the catalog holds in each `(period, repeat count)` stratum.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StratumCounts {
    counts: HashMap<(u8, u64), u64>,
}

impl StratumCounts {
    /// Count one locus.
    pub fn count(&mut self, period: u8, repeat_count: u64) {
        *self.counts.entry((period, repeat_count)).or_insert(0) += 1;
    }

    /// How many loci that stratum holds.
    pub fn get(&self, period: u8, repeat_count: u64) -> u64 {
        self.counts
            .get(&(period, repeat_count))
            .copied()
            .unwrap_or(0)
    }

    /// Every stratum with at least one locus, in `(period, repeat count)` order — sorted so
    /// that a report or a test reads the same way on every machine.
    pub fn iter_sorted(&self) -> Vec<((u8, u64), u64)> {
        let mut out: Vec<((u8, u64), u64)> = self.counts.iter().map(|(k, v)| (*k, *v)).collect();
        out.sort_unstable();
        out
    }

    /// Loci counted, all strata.
    pub fn total(&self) -> u64 {
        self.counts.values().sum()
    }

    /// The counts a census file wrote down, read back.
    ///
    /// **The inverse of [`iter_sorted`](Self::iter_sorted) and nothing else.** Counting is what
    /// builds these; this exists so a set of counts that was counted can survive a round trip
    /// through bytes without being recounted from loci the file does not carry.
    pub fn from_counted(entries: impl IntoIterator<Item = ((u8, u64), u64)>) -> Self {
        Self {
            counts: entries.into_iter().collect(),
        }
    }
}

/// Up to `cap` loci from each stratum: the ones whose hash is lowest.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StratumSample {
    kept: HashMap<(u8, u64), Vec<SsrSegment>>,
}

impl StratumSample {
    /// The loci kept for that stratum, in coordinate order.
    pub fn get(&self, period: u8, repeat_count: u64) -> &[SsrSegment] {
        self.kept
            .get(&(period, repeat_count))
            .map(Vec::as_slice)
            .unwrap_or(&[])
    }

    /// Every stratum sampled, in `(period, repeat count)` order.
    pub fn iter_sorted(&self) -> Vec<((u8, u64), &[SsrSegment])> {
        let mut keys: Vec<(u8, u64)> = self.kept.keys().copied().collect();
        keys.sort_unstable();
        keys.into_iter()
            .map(|k| (k, self.kept[&k].as_slice()))
            .collect()
    }

    /// Loci kept, all strata.
    pub fn total(&self) -> usize {
        self.kept.values().map(Vec::len).sum()
    }
}

/// Accumulates the `cap` lowest-hash loci per stratum, and the full counts beside them.
///
/// **One pass, bounded state, order-independent.** Each locus is hashed once and offered to
/// its stratum's heap, which never holds more than `cap`; nothing is sorted and no locus is
/// visited twice. The result does not depend on the order loci arrive in, which is what lets
/// a sharded enumeration merge by taking the lowest `cap` of the union
/// (`parameter_prepass_joint_loci.md` §3.2).
pub struct StratumSampler {
    cap: usize,
    seed: u64,
    counts: StratumCounts,
    /// Per stratum, the kept loci with their hashes. Kept as a plain vector rather than a
    /// `BinaryHeap` because `cap` is small and the max is scanned only when the vector is
    /// full — the cost is `cap` comparisons on the rare replacement, against a heap's
    /// allocation per stratum.
    kept: HashMap<(u8, u64), Vec<(u64, SsrSegment)>>,
}

impl StratumSampler {
    /// A sampler keeping `cap` loci per stratum, hashing with `seed`.
    ///
    /// The seed is the caller's: two runs at the same seed must keep the identical loci, and
    /// that is a property of the run rather than of the file.
    pub fn new(cap: usize, seed: u64) -> Self {
        Self {
            cap,
            seed,
            counts: StratumCounts::default(),
            kept: HashMap::new(),
        }
    }

    /// Offer one locus, of the given stratum.
    pub fn offer(&mut self, period: u8, repeat_count: u64, locus: SsrSegment) {
        self.counts.count(period, repeat_count);
        if self.cap == 0 {
            return;
        }

        let value = hash_locus(locus.chrom(), locus.start(), self.seed);
        let bucket = self.kept.entry((period, repeat_count)).or_default();
        if bucket.len() < self.cap {
            bucket.push((value, locus));
            return;
        }

        // Full: replace the highest hash, but only if this one is lower. Ties keep the
        // incumbent, which is what makes the result order-independent — two loci with the
        // same hash are the same locus, since the hash is over its coordinates.
        let (worst_index, worst_value) = bucket
            .iter()
            .enumerate()
            .map(|(i, (v, _))| (i, *v))
            .max_by_key(|(_, v)| *v)
            .expect("a full bucket is not empty");
        if value < worst_value {
            bucket[worst_index] = (value, locus);
        }
    }

    /// The counts and the sample, the counts covering every locus offered and the sample the
    /// `cap` lowest hashes of each stratum.
    pub fn finish(self) -> (StratumCounts, StratumSample) {
        let kept = self
            .kept
            .into_iter()
            .map(|(stratum, mut loci)| {
                // Coordinate order, so a caller reads loci the way it reads everything else.
                loci.sort_by_key(|(_, locus)| locus.start());
                (stratum, loci.into_iter().map(|(_, locus)| locus).collect())
            })
            .collect();
        (self.counts, StratumSample { kept })
    }
}

/// The value a locus is selected by: a hash of **where it is**, not of what it holds.
///
/// Position, not content, so that two loci with the same motif and length are independent
/// draws; and a fixed hash rather than a random number generator, so the same reference,
/// seed and cap keep the identical loci on every machine and in every run.
fn hash_locus(chrom: &str, start: u64, seed: u64) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut hasher = xxhash_rust::xxh3::Xxh3::with_seed(seed);
    chrom.hash(&mut hasher);
    start.hash(&mut hasher);
    hasher.finish()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::types::Motif;

    fn locus(chrom: &str, start: u64, copies: u64, period: u8) -> SsrSegment {
        let motif = match period {
            1 => b"A".to_vec(),
            2 => b"AT".to_vec(),
            _ => b"CAG".to_vec(),
        };
        SsrSegment::new(
            chrom.to_string().into_boxed_str(),
            start,
            start + copies * u64::from(period) - 1,
            Motif::new(&motif).expect("valid"),
            1.0,
        )
        .expect("valid locus")
    }

    /// A stratum smaller than the cap keeps everything it has.
    #[test]
    fn a_small_stratum_keeps_every_locus() {
        let mut sampler = StratumSampler::new(10, 7);
        for i in 0..4 {
            sampler.offer(3, 8, locus("chr1", 100 + i * 100, 8, 3));
        }
        let (counts, sample) = sampler.finish();
        assert_eq!(counts.get(3, 8), 4);
        assert_eq!(sample.get(3, 8).len(), 4);
    }

    /// A stratum larger than the cap keeps exactly `cap`, and the counts still see them all
    /// — the pre-pass reweights by the total, so a truncated count would bias it.
    #[test]
    fn a_large_stratum_keeps_exactly_the_cap_and_counts_them_all() {
        let mut sampler = StratumSampler::new(5, 7);
        for i in 0..100 {
            sampler.offer(2, 12, locus("chr1", 1_000 + i * 50, 12, 2));
        }
        let (counts, sample) = sampler.finish();
        assert_eq!(counts.get(2, 12), 100, "every locus is counted");
        assert_eq!(sample.get(2, 12).len(), 5, "exactly the cap is kept");
    }

    /// The same seed keeps the identical loci; a different seed keeps a different set. The
    /// first is what makes a run reproducible, the second that the sample is not the first
    /// `cap` in disguise.
    #[test]
    fn the_sample_is_stable_under_a_seed_and_moves_with_it() {
        let starts: Vec<u64> = (0..60).map(|i| 1_000 + i * 37).collect();
        let sample_at = |seed: u64| {
            let mut sampler = StratumSampler::new(6, seed);
            for &start in &starts {
                sampler.offer(3, 9, locus("chr1", start, 9, 3));
            }
            let (_, sample) = sampler.finish();
            sample
                .get(3, 9)
                .iter()
                .map(|l| l.start())
                .collect::<Vec<_>>()
        };

        assert_eq!(sample_at(11), sample_at(11), "same seed, same loci");
        assert_ne!(
            sample_at(11),
            sample_at(12),
            "a different seed draws a different sample"
        );
    }

    /// **Order-independence**, which is what lets a sharded enumeration merge: the same loci
    /// offered in a different order give the same sample, and merging two halves by taking
    /// the lowest `cap` of the union gives it too.
    #[test]
    fn the_sample_does_not_depend_on_the_order_loci_arrive_in() {
        let starts: Vec<u64> = (0..80).map(|i| 500 + i * 23).collect();
        let sampled = |order: &[u64]| {
            let mut sampler = StratumSampler::new(7, 99);
            for &start in order {
                sampler.offer(2, 11, locus("chr2", start, 11, 2));
            }
            let (_, sample) = sampler.finish();
            sample
                .get(2, 11)
                .iter()
                .map(|l| l.start())
                .collect::<Vec<_>>()
        };

        let forwards = sampled(&starts);
        let mut backwards_order = starts.clone();
        backwards_order.reverse();
        assert_eq!(forwards, sampled(&backwards_order), "reversed order");

        // Two shards, merged by re-offering each half's survivors: the union's lowest `cap`.
        let (left, right) = starts.split_at(starts.len() / 2);
        let mut merged = StratumSampler::new(7, 99);
        for start in sampled(left).into_iter().chain(sampled(right)) {
            merged.offer(2, 11, locus("chr2", start, 11, 2));
        }
        let (_, merged_sample) = merged.finish();
        let merged_starts: Vec<u64> = merged_sample.get(2, 11).iter().map(|l| l.start()).collect();
        assert_eq!(
            merged_starts, forwards,
            "merging two shards is taking the lowest cap of the union"
        );
    }

    /// The sample is drawn per stratum, not across them: a rare stratum is not crowded out
    /// by a common one, which is the whole reason the pre-pass samples this way.
    #[test]
    fn each_stratum_is_sampled_on_its_own() {
        let mut sampler = StratumSampler::new(3, 5);
        for i in 0..50 {
            sampler.offer(2, 10, locus("chr1", 1_000 + i * 40, 10, 2));
        }
        for i in 0..4 {
            sampler.offer(3, 20, locus("chr1", 50_000 + i * 100, 20, 3));
        }
        let (counts, sample) = sampler.finish();

        assert_eq!(counts.get(2, 10), 50);
        assert_eq!(counts.get(3, 20), 4);
        assert_eq!(sample.get(2, 10).len(), 3, "the common stratum is capped");
        assert_eq!(
            sample.get(3, 20).len(),
            3,
            "the rare one is capped the same"
        );
    }
}
