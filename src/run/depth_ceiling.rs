//! **The depth ceiling: a locus where any sample has a read group deeper than this is not
//! called**, and none of its evidence is held (`--max-read-group-depth`).
//!
//! # What it is for
//!
//! A position where thousands of reads pile up in a sample that averages three is almost always
//! a collapsed repeat: many copies of a sequence in the genome, mapped onto one copy in the
//! reference. Nothing called there can be trusted, and holding the evidence costs memory in
//! proportion to depth × records × samples — which is what stopped the 2,169-sample tomato run
//! at SL4.0ch00:1,504,001, at about 48 MB a sample after the stored read lists were already
//! reduced to their changes (`doc/devel/reports/psp_pileup_memory_2026-09-28/report.md`).
//!
//! So a sample's record over the ceiling is not held. Its source hands the merge only where it
//! lay ([`LocusSummary::over_depth_ceiling`](super::cohort_merge::LocusSummary::over_depth_ceiling)),
//! and **the merge builds no cohort locus over that ground for any sample** (owner, 2026-09-29):
//! a locus that one sample's reads say is a pile-up is not one the others can call either.
//!
//! **This is a different thing from the read cap** (`--max-reads-per-position`). The cap keeps
//! a locus and uses at most that many of a read group's reads there; the ceiling drops the locus.
//!
//! # How deep a read group is, from what a record stores
//!
//! A record stores its reads as observations, each tagged with its read group and a count, plus
//! two counts that belong to no read group: the reads that covered the locus but showed nothing,
//! and the reads the walk's own cap discarded. **Both are counted with the record's deepest read
//! group.** The discarded ones come from a read group the walk's cap was cutting, which is the
//! deepest one; the silent ones cannot be placed, and adding them where they tip the answer
//! soonest errs toward dropping. For a sample with one read group, which is the usual shape,
//! this is simply the record's depth.
//!
//! So a record is over the ceiling when its deepest read group's observed reads, plus those two
//! counts, exceed it.

use crate::locus_generation::SampleLocusObservations;

/// The default of `--max-read-group-depth`: 1,000 reads of one read group at one locus.
///
/// **The owner's choice (2026-09-29).** It is far above any honest depth at the coverages the
/// caller is for — a 30× sample is thirty times below it, a 300× sample three times — and far
/// below the pile-ups it exists for, which ran to a median of 5,669 reads a position in 3×
/// samples.
pub const DEFAULT_MAX_READ_GROUP_DEPTH: u32 = 1_000;

/// The most reads one read group may have at a locus before the locus is dropped for the whole
/// cohort. See the module documentation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DepthCeiling {
    reads: u32,
}

impl Default for DepthCeiling {
    fn default() -> Self {
        Self::at(DEFAULT_MAX_READ_GROUP_DEPTH)
    }
}

impl DepthCeiling {
    /// A ceiling of `reads` reads of one read group.
    #[must_use]
    pub fn at(reads: u32) -> Self {
        Self { reads }
    }

    /// The most reads a read group may have.
    #[must_use]
    pub fn reads(self) -> u32 {
        self.reads
    }

    /// Whether `record` has a read group deeper than the ceiling.
    #[must_use]
    pub fn is_exceeded_by(self, record: &SampleLocusObservations) -> bool {
        let unplaced =
            u64::from(record.reads_without_observation) + u64::from(record.reads_discarded_by_cap);
        let ceiling = u64::from(self.reads);
        if unplaced > ceiling {
            return true;
        }
        // **The ordinary record leaves here, having allocated nothing**: when every read it holds,
        // summed over its read groups, is within the ceiling, no one read group can be over it.
        let observed: u64 = record
            .observations
            .iter()
            .map(|observation| u64::from(observation.num_obs))
            .sum();
        if observed + unplaced <= ceiling {
            return false;
        }
        deepest_read_group(record) + unplaced > ceiling
    }

    /// **Whether a stored record could be over the ceiling, from what its head walk already
    /// knows** — so that the reader decodes only the records that could, and not the ninety-nine
    /// in a hundred it otherwise never decodes.
    ///
    /// `spanning_fragments` is how many read fragments the file lists as covering the record's
    /// position (a pair of mates is one fragment) and `discarded` the reads the walk's cap
    /// discarded there. A fragment is at most two reads, so twice the fragments plus the
    /// discarded reads bounds the reads a record can hold.
    ///
    /// **Except where the file lists no fragments at all**, which is a repeat tract: its
    /// observations are counted but name no reads, so nothing bounds them from the head, and
    /// such a record is always decoded. Measured over ch00:1–1,600,000 of one median-sized
    /// tomato psp (SRS11350780, about 3×, both of that stretch's pile-ups): of 1,174,134
    /// records, the 334 with no fragments listed were exactly its repeat tracts; 175 of them held
    /// more reads than twice their fragments, 14 more than 100. Every one of the other 1,173,800
    /// held no more reads than the bound.
    #[must_use]
    pub fn may_be_exceeded_at(self, spanning_fragments: usize, discarded: u32) -> bool {
        spanning_fragments == 0
            || 2 * spanning_fragments as u64 + u64::from(discarded) > u64::from(self.reads)
    }
}

/// The observed reads of `record`'s deepest read group.
///
/// **A short list scanned per group rather than a map**: a sample has a handful of read groups,
/// and this runs only on the records deep enough to need it.
fn deepest_read_group(record: &SampleLocusObservations) -> u64 {
    let mut per_group: Vec<(u32, u64)> = Vec::new();
    for observation in &record.observations {
        let group = observation.read_group.get();
        let reads = u64::from(observation.num_obs);
        match per_group.iter_mut().find(|(seen, _)| *seen == group) {
            Some((_, total)) => *total += reads,
            None => per_group.push((group, reads)),
        }
    }
    per_group.iter().map(|(_, reads)| *reads).max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::locus_generation::{LocusKind, ReadWitness, SequenceObservation};
    use crate::types::{ContigId, GenomeRegion, Position, ReadGroupId, SummedLogError};

    fn observation(group: u32, reads: u32) -> SequenceObservation {
        SequenceObservation {
            bases: b"A".to_vec().into_boxed_slice(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(group),
            num_obs: reads,
            num_fwd: 0,
            q_sum: SummedLogError::from_nats(0.0),
            mapq_sum: 0,
            mapq_sum_sq: 0,
            placed_left: 0,
            chain_ids: Vec::new(),
        }
    }

    fn record(
        observations: Vec<SequenceObservation>,
        silent: u32,
        discarded: u32,
    ) -> SampleLocusObservations {
        SampleLocusObservations {
            region: GenomeRegion {
                contig: ContigId(0),
                start: Position(10),
                end: Position(10),
            },
            reference_bases: b"A".to_vec().into_boxed_slice(),
            observations,
            reads_without_observation: silent,
            reads_discarded_by_cap: discarded,
            kind: LocusKind::Generic,
        }
    }

    #[test]
    fn one_read_group_is_judged_on_the_records_whole_depth() {
        let ceiling = DepthCeiling::at(100);
        assert!(!ceiling.is_exceeded_by(&record(
            vec![observation(0, 60), observation(0, 40)],
            0,
            0
        )));
        assert!(ceiling.is_exceeded_by(&record(
            vec![observation(0, 60), observation(0, 41)],
            0,
            0
        )));
    }

    #[test]
    fn read_groups_are_judged_apart() {
        // 150 reads in all, but no group has more than 80.
        let ceiling = DepthCeiling::at(100);
        let spread = record(vec![observation(0, 70), observation(1, 80)], 0, 0);
        assert!(!ceiling.is_exceeded_by(&spread));
        let one_deep = record(vec![observation(0, 30), observation(1, 101)], 0, 0);
        assert!(ceiling.is_exceeded_by(&one_deep));
    }

    #[test]
    fn discarded_and_silent_reads_count_with_the_deepest_group() {
        let ceiling = DepthCeiling::at(100);
        // The walk kept 90 of group 1 and discarded 20: group 1 had 110.
        let capped = record(vec![observation(0, 30), observation(1, 90)], 0, 20);
        assert!(ceiling.is_exceeded_by(&capped));
        let silent = record(vec![observation(0, 30), observation(1, 90)], 11, 0);
        assert!(ceiling.is_exceeded_by(&silent));
        let within = record(vec![observation(0, 30), observation(1, 90)], 5, 5);
        assert!(!ceiling.is_exceeded_by(&within));
    }

    /// Fifty fragments are at most a hundred reads; fifty-one could be more. A record listing no
    /// fragments is a repeat tract, which nothing in the head bounds.
    #[test]
    fn the_head_rules_out_only_records_it_can_bound() {
        let ceiling = DepthCeiling::at(100);
        assert!(!ceiling.may_be_exceeded_at(50, 0));
        assert!(ceiling.may_be_exceeded_at(51, 0));
        assert!(ceiling.may_be_exceeded_at(1, 99));
        assert!(ceiling.may_be_exceeded_at(0, 0));
    }
}
