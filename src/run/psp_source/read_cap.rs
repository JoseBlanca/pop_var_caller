//! **The read cap applied again as a stored record is decoded**, for psps written under a
//! looser one.
//!
//! `--max-reads-per-position` caps the reads of one read group at one position when the
//! observations are built ([`WalkerConfig`](crate::locus_generation::pileup)). Psps written
//! before 2026-09-29 were built under 8,000 reads a *sample*, and regenerating the 2,169-sample
//! tomato cohort (8.6 TB) to lower that is not an option, so a run reading them caps them here
//! instead.
//!
//! # Why it thins rather than drops reads
//!
//! **A psp does not store reads one by one.** An observation is one allele seen in one read
//! group, with the count of reads that showed it, how many were on the forward strand, how
//! many started left of the locus, and the *sums* of their error and mapping quality
//! ([`SequenceObservation`]). The sums cannot be split back into reads, so the cap cannot be
//! re-run as the walk runs it. What can be done is to keep the same share of every
//! observation in the read group: the counts and the sums shrink by that share, so each
//! allele's mean error and mean mapping quality are unchanged and the balance between alleles
//! holds to within rounding (owner, 2026-09-29).
//!
//! # Which reads, where the record names them
//!
//! **The share is taken as whole reads**, so the same reads are kept at neighbouring positions:
//! a read is kept when the hash of its identifier falls in the lowest `share` of the hash's
//! range. A fixed cutoff rather than "the lowest-hashing reads, as many as the share allows",
//! because the latter moves with each position's read count and would keep a read at one
//! position and drop it at the next; with the cutoff a read is kept wherever the share is the
//! same, and a read kept at a smaller share is kept at every larger one. The price is that the
//! count kept is the cap only on average — at 1,000 it scatters by about 30 either way. Each
//! observation's count and sums are then scaled by the fraction of its own reads that were
//! kept. **Repeat loci name no reads** ([`SequenceObservation::chain_ids`] is empty there), and
//! for them every observation keeps the group's share of its count directly.
//!
//! **It is not the walk's cap.** The walk keeps reads by a hash of their *name*, which a psp
//! does not store; this keeps them by a hash of their identifier within the file. Calling from
//! thinned psps therefore does not reproduce calling from the alignments at the same cap — it
//! is close, not equal. A psp written under a cap no looser than the run's is never thinned,
//! so fresh files lose nothing (see [`needs_thinning`]).

use crate::locus_generation::{SampleLocusObservations, SequenceObservation};
use crate::psp::{Header, ParameterValue};
use crate::run::gatherer::MAX_READS_PER_POSITION_KEY;
use crate::types::SummedLogError;

/// Whether a run capping at `run_cap` reads per read group has to thin the records of the psp
/// whose header is `header`.
///
/// **Only when the file was written under a looser cap.** A file that records one no looser
/// than the run's already holds at most that many reads of a read group at a position, so
/// thinning it would change nothing but the reads it names at a multi-position record — and
/// would make calling from fresh psps differ from calling from the alignments. A file that
/// records none predates the key and was written under 8,000 reads a sample, so it is thinned.
#[must_use]
pub fn needs_thinning(header: &Header, run_cap: u32) -> bool {
    match header.writer.parameters.get(MAX_READS_PER_POSITION_KEY) {
        Some(ParameterValue::Integer(file_cap)) => *file_cap > i64::from(run_cap),
        _ => true,
    }
}

/// What thinning one record removed.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Thinned {
    /// Reads removed from the record's observations, summed over its read groups.
    pub reads: u64,
}

/// Thin every read group of `record` that has more than `cap` reads down to about `cap`,
/// adding the reads removed to its `reads_discarded_by_cap`. See the module documentation.
pub fn thin_to_read_cap(record: &mut SampleLocusObservations, cap: u32) -> Thinned {
    // **The ordinary record leaves here, having allocated nothing.** This runs on every record a
    // run decodes from an older psp, and almost none of them hold more reads than the cap even
    // summed over every read group — in which case no single group can.
    let reads: u64 = record
        .observations
        .iter()
        .map(|o| u64::from(o.num_obs))
        .sum();
    if reads <= u64::from(cap) {
        return Thinned::default();
    }
    let mut groups: Vec<_> = record.observations.iter().map(|o| o.read_group).collect();
    groups.sort_unstable();
    groups.dedup();
    let mut removed = 0u64;
    for group in groups {
        let reads: u64 = record
            .observations
            .iter()
            .filter(|o| o.read_group == group)
            .map(|o| u64::from(o.num_obs))
            .sum();
        if reads <= u64::from(cap) {
            continue;
        }
        let share = f64::from(cap) / reads as f64;
        let cutoff = (share * u64::MAX as f64) as u64;
        for observation in record
            .observations
            .iter_mut()
            .filter(|o| o.read_group == group)
        {
            let before = observation.num_obs;
            let fraction = if observation.chain_ids.is_empty() {
                share
            } else {
                let named = observation.chain_ids.len();
                observation.chain_ids.retain(|id| splitmix64(*id) < cutoff);
                observation.chain_ids.len() as f64 / named as f64
            };
            scale(observation, fraction);
            removed += u64::from(before - observation.num_obs);
        }
    }
    record.observations.retain(|o| o.num_obs > 0);
    record.reads_discarded_by_cap = record
        .reads_discarded_by_cap
        .saturating_add(u32::try_from(removed).unwrap_or(u32::MAX));
    Thinned { reads: removed }
}

/// Keep `fraction` of `observation`'s reads: every count and sum scaled by it, the counts
/// rounded and the parts never above the whole.
fn scale(observation: &mut SequenceObservation, fraction: f64) {
    let count = |value: u32| (f64::from(value) * fraction).round() as u32;
    let kept = count(observation.num_obs);
    if observation.num_obs == 0 {
        return;
    }
    // Each sum is scaled by the fraction of reads actually kept, after rounding, so its mean
    // over the kept reads is the mean it had over all of them.
    let exact = f64::from(kept) / f64::from(observation.num_obs);
    observation.num_fwd = count(observation.num_fwd).min(kept);
    observation.placed_left = count(observation.placed_left).min(kept);
    observation.q_sum = SummedLogError::from_nats(observation.q_sum.nats() * exact);
    observation.mapq_sum = (f64::from(observation.mapq_sum) * exact).round() as u32;
    observation.mapq_sum_sq = (observation.mapq_sum_sq as f64 * exact).round() as u64;
    observation.num_obs = kept;
}

/// The splitmix64 finaliser: a bijection that spreads consecutive identifiers across the whole
/// range, so "lowest hash" is not "lowest identifier" — which would keep the reads that arrived
/// first, a sample tilted towards early alignment starts.
fn splitmix64(value: u64) -> u64 {
    let mut z = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::locus_generation::ReadWitness;
    use crate::types::ReadGroupId;

    fn observation(group: u32, bases: &[u8], ids: std::ops::Range<u64>) -> SequenceObservation {
        let reads = (ids.end - ids.start) as u32;
        SequenceObservation {
            bases: bases.into(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(group),
            num_obs: reads,
            num_fwd: reads / 2,
            q_sum: SummedLogError::from_nats(-0.01 * f64::from(reads)),
            mapq_sum: 60 * reads,
            mapq_sum_sq: 3_600 * u64::from(reads),
            placed_left: reads / 4,
            chain_ids: ids.collect(),
        }
    }

    fn record(observations: Vec<SequenceObservation>) -> SampleLocusObservations {
        SampleLocusObservations {
            region: crate::types::GenomeRegion {
                contig: crate::types::ContigId(0),
                start: crate::types::Position(1),
                end: crate::types::Position(1),
            },
            reference_bases: b"A".to_vec().into_boxed_slice(),
            observations,
            reads_without_observation: 0,
            reads_discarded_by_cap: 0,
            kind: crate::locus_generation::LocusKind::Generic,
        }
    }

    /// **A group over the cap keeps about the cap, each allele its share, and the means hold.**
    /// 4,000 reads in one group, 3,000 showing A and 1,000 showing C, capped at 1,000.
    #[test]
    fn a_group_over_the_cap_keeps_each_alleles_share_and_its_means() {
        let mut thinned = record(vec![
            observation(0, b"A", 0..3_000),
            observation(0, b"C", 3_000..4_000),
        ]);
        let removed = thin_to_read_cap(&mut thinned, 1_000);
        let total: u32 = thinned.observations.iter().map(|o| o.num_obs).sum();
        assert!((900..=1_100).contains(&total), "kept {total}");
        assert_eq!(removed.reads, u64::from(4_000 - total));
        assert_eq!(thinned.reads_discarded_by_cap, 4_000 - total);
        let a = &thinned.observations[0];
        let share_of_a = f64::from(a.num_obs) / f64::from(total);
        assert!(
            (share_of_a - 0.75).abs() < 0.03,
            "A is {share_of_a} of the kept reads"
        );
        for kept in &thinned.observations {
            assert_eq!(kept.chain_ids.len() as u32, kept.num_obs);
            assert_eq!(
                kept.mapq_sum,
                60 * kept.num_obs,
                "the mean mapping quality holds"
            );
            let mean_error = kept.q_sum.nats() / f64::from(kept.num_obs);
            assert!(
                (mean_error + 0.01).abs() < 1e-3,
                "the mean error holds: {mean_error}"
            );
            assert!(kept.num_fwd <= kept.num_obs && kept.placed_left <= kept.num_obs);
        }
    }

    /// **Neighbouring positions keep the same reads**, which is what lets a later step chain a
    /// read's alleles across them: two positions with the same share keep exactly the same
    /// reads of those they share.
    #[test]
    fn the_same_reads_are_kept_at_neighbouring_positions() {
        let mut here = record(vec![observation(0, b"A", 0..4_000)]);
        let mut next = record(vec![observation(0, b"A", 10..4_010)]);
        thin_to_read_cap(&mut here, 1_000);
        thin_to_read_cap(&mut next, 1_000);
        let here: std::collections::HashSet<_> = here.observations[0].chain_ids.iter().collect();
        let shared_kept = next.observations[0]
            .chain_ids
            .iter()
            .filter(|id| **id < 4_000 && here.contains(id))
            .count();
        let next_kept_among_shared = next.observations[0]
            .chain_ids
            .iter()
            .filter(|id| **id < 4_000)
            .count();
        assert_eq!(
            shared_kept, next_kept_among_shared,
            "a read kept at the second position is kept at the first"
        );
    }

    /// **Each read group is capped on its own**, and a group within the cap is untouched.
    #[test]
    fn each_read_group_is_capped_on_its_own() {
        let mut thinned = record(vec![
            observation(0, b"A", 0..3_000),
            observation(1, b"A", 3_000..3_500),
        ]);
        thin_to_read_cap(&mut thinned, 1_000);
        let of = |group: u32| -> u32 {
            thinned
                .observations
                .iter()
                .filter(|o| o.read_group == ReadGroupId(group))
                .map(|o| o.num_obs)
                .sum()
        };
        assert!((900..=1_100).contains(&of(0)));
        assert_eq!(of(1), 500, "the group within the cap keeps every read");
    }

    /// **A repeat locus names no reads**, so its observations keep the group's share directly.
    #[test]
    fn observations_that_name_no_reads_keep_the_groups_share() {
        let mut unnamed = observation(0, b"ATATAT", 0..3_000);
        unnamed.chain_ids.clear();
        let mut other = observation(0, b"ATATATAT", 3_000..4_000);
        other.chain_ids.clear();
        let mut thinned = record(vec![unnamed, other]);
        thin_to_read_cap(&mut thinned, 1_000);
        assert_eq!(thinned.observations[0].num_obs, 750);
        assert_eq!(thinned.observations[1].num_obs, 250);
    }

    /// **Nothing moves under the cap.**
    #[test]
    fn a_record_within_the_cap_is_untouched() {
        let original = record(vec![
            observation(0, b"A", 0..600),
            observation(0, b"C", 600..1_000),
        ]);
        let mut thinned = original.clone();
        assert_eq!(thin_to_read_cap(&mut thinned, 1_000), Thinned::default());
        assert_eq!(thinned, original);
    }
}
