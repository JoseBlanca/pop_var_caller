//! What one builder's walk costs as the cohort grows — the measurement behind the merge's
//! k-way merge.
//!
//! The walk takes observations in genome order, and to do that it must answer, for every
//! observation, "which sample comes next?". This times one 20-base building region with every
//! sample carrying a record at every position of it — which is what the generic mint produces
//! wherever a sample has reads — at cohort sizes spanning the range the caller commits to.
//!
//! **It reports the median of several repeats and their spread**, because the machine's own
//! swing is large enough to be mistaken for a code change: repeated runs of one unchanged
//! binary at 3,000 samples gave 3,342 µs, 3,412 µs and 4,438 µs. One mean cannot tell that
//! from a third off.
//!
//! Run in release: `cargo run --release --example ng_cohort_merge_walk_cost`.

use std::time::Instant;

use pop_var_caller::ng::locus_generation::{
    LocusKind, ReadWitness, SampleLocusObservations, SequenceObservation,
};
use pop_var_caller::ng::run::cohort_merge::build::build_region;
use pop_var_caller::ng::run::cohort_merge::{MaxCohortLocusSpan, MinAltReads};
use pop_var_caller::ng::types::{ContigId, GenomeRegion, Position, ReadGroupId};

/// One sample's record at one position, showing the reference to three reads.
fn record_at(position: u64, observed: &[u8]) -> SampleLocusObservations {
    SampleLocusObservations {
        region: GenomeRegion {
            contig: ContigId(0),
            start: Position(position),
            end: Position(position),
        },
        reference_bases: b"A".to_vec(),
        observations: vec![SequenceObservation {
            bases: observed.to_vec(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: 3,
            num_fwd: 3,
            q_sum: pop_var_caller::ng::types::SummedLogError::from_nats(-6.0),
            mapq_sum: 180,
            mapq_sum_sq: 10_800,
            placed_left: 0,
            chain_ids: vec![1, 2, 3],
        }],
        reads_without_observation: 0,
        reads_discarded_by_cap: 0,
        kind: LocusKind::Generic,
    }
}

fn main() {
    let region = GenomeRegion {
        contig: ContigId(0),
        start: Position(1),
        end: Position(20),
    };
    // One column, not two: a region where one sample carries a substitution measured within
    // the machine's own swing of an all-reference one at every cohort size, so quoting both
    // would imply a difference the measurement does not carry.
    println!("cohort, median_us, min_us, max_us");

    for samples in [1usize, 10, 63, 250, 1000, 3000] {
        // Every sample covers every position of the region, all agreeing with the reference.
        let quiet: Vec<Vec<SampleLocusObservations>> = (0..samples)
            .map(|_| (1..=20).map(|at| record_at(at, b"A")).collect())
            .collect();
        let time = |cohort: &[Vec<SampleLocusObservations>]| {
            let slices: Vec<&[SampleLocusObservations]> =
                cohort.iter().map(Vec::as_slice).collect();
            // Warm, then seven repeats of enough walks that each timed block is long
            // beside the clock and the machine's own jitter.
            let _ = build_region(
                region,
                &slices,
                MaxCohortLocusSpan::DEFAULT,
                MinAltReads::DEFAULT,
            );
            let walks = if samples >= 1000 { 30 } else { 300 };
            let mut each_repeat: Vec<f64> = (0..7)
                .map(|_| {
                    let started = Instant::now();
                    for _ in 0..walks {
                        let outcome = build_region(
                            region,
                            &slices,
                            MaxCohortLocusSpan::DEFAULT,
                            MinAltReads::DEFAULT,
                        );
                        std::hint::black_box(&outcome);
                    }
                    started.elapsed().as_secs_f64() * 1e6 / f64::from(walks)
                })
                .collect();
            each_repeat.sort_by(f64::total_cmp);
            (
                each_repeat[each_repeat.len() / 2],
                each_repeat[0],
                each_repeat[each_repeat.len() - 1],
            )
        };

        let (median, fastest, slowest) = time(&quiet);
        println!("{samples}, {median:.2}, {fastest:.2}, {slowest:.2}");
    }
}
