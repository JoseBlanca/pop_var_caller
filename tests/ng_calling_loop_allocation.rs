//! **Nothing allocates inside the calling loop's passes** — `doc/devel/ng/spec/calling_em_loop.md`
//! §13's test 7, counted rather than inferred.
//!
//! # Why this is a test binary of its own
//!
//! A `#[global_allocator]` counts **the whole process**, and `cargo test` runs a binary's tests
//! in parallel by default — so a counter read inside one test of the library suite would include
//! whatever the tests running beside it allocated. Here the binary holds one test, so the only
//! allocations between the two readings are the ones this test caused.
//!
//! # Why it needs a feature, and why it costs the crate nothing
//!
//! The allocator is `dhat`'s, which is already this repository's heap-profiling dependency
//! (`Cargo.toml`'s `dhat-heap`). **No `unsafe` is written here or anywhere in this crate**:
//! `#[global_allocator]` is a safe attribute and the `unsafe impl GlobalAlloc` behind it is
//! `dhat`'s own, so `src/lib.rs`'s `#![forbid(unsafe_code)]` stands untouched. Without the
//! feature this file compiles to nothing, and the library suite's
//! `no_buffer_of_the_scratch_moves_or_grows_however_many_passes_the_loop_takes` is the cheap
//! guard that runs on every build.
//!
//! ```text
//! ./scripts/dev.sh cargo test --test ng_calling_loop_allocation --features dhat-heap
//! ```
//!
//! # What it would catch that the library suite cannot
//!
//! A temporary — a `Vec` allocated and dropped inside a pass — leaves no trace in any scratch
//! buffer, so the pointer-and-length fingerprint the library suite compares cannot see it.
//! Measured on the uncontaminated arm: its two runs allocate **8 blocks each**, and one
//! `Vec::with_capacity` added to the frequency loop's seeded pass takes them apart — **8 against
//! 10**, one block per extra pass — while every one of the library's 4,694 tests stayed green.
//! Measured again on the contaminated arm (2026-08-26): its two runs allocate **8 blocks each**,
//! and one `Vec::with_capacity` added beside the per-pass table assembly gives **8 against 11**,
//! again one block per extra pass.

#![cfg(feature = "dhat-heap")]

use pop_var_caller::ng::calling::genotype_prior::{
    MarginalizedDirichletPrior, SeedRegime, SpectrumSeed,
};
use pop_var_caller::ng::calling::inference::summarise_condition::SummariseConditionLoop;
use pop_var_caller::ng::calling::inference::{CallingLoopConfig, LocusGenotyper};
use pop_var_caller::ng::calling::likelihood::ssr_emission::{
    StutterSubstitutionEmission, StutterSubstitutionScratch,
};
use pop_var_caller::ng::calling::{
    CallingScratch, CandidateAlleles, ContaminationView, FrozenParameters, GenericLocusSample,
    GenericObservation, GenericSampleEvidence, LocusEvidence, ReadGroupCalibration,
};
use pop_var_caller::ng::locus_generation::LocusKind;
use pop_var_caller::ng::parameter_estimation::joint::contamination::ContaminationSource;
use pop_var_caller::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
use pop_var_caller::ng::parameter_estimation::joint::stratum_fits::StratumFits;
use pop_var_caller::ng::types::{
    AlleleId, ContigId, GenomeRegion, InbreedingF, Ploidy, Position, ReadGroupId,
};
use std::num::NonZeroU32;

#[global_allocator]
static COUNTING_ALLOCATOR: dhat::Alloc = dhat::Alloc;

/// One `(allele, read group)` row of a sample's evidence.
fn observation(allele: u16, num_reads: u32) -> GenericObservation {
    GenericObservation {
        allele: AlleleId(allele),
        read_group: ReadGroupId(0),
        num_reads,
        q_sum: -3.0 * f64::from(num_reads),
        forward_reads: num_reads / 2,
        placed_left_reads: num_reads / 2,
    }
}

/// **The same locus at two pass counts allocates the same number of blocks — with and without
/// contamination.**
///
/// Three samples over three candidate alleles, called once with the pass cap at two and once at
/// the shipped default. The loop's buffers all belong to the worker's scratch and are sized once
/// per locus, so the extra passes must cost nothing — and "nothing" here is a count, not an
/// inference from a pointer.
///
/// **The contaminated arm is the one that does per-pass work inside the loop.** With no fraction
/// fitted the loop reads a table nothing touches; with one fitted it refills the batch copies,
/// refills one sample's contaminant frequencies per row and assembles the whole table again, at
/// every pass. Those three buffers are sized outside the loop and written inside it, which is
/// exactly the shape a temporary hides in.
///
/// **Both arms live in this one test rather than in two**, because a `#[global_allocator]` counts
/// the whole process and `cargo test` runs a binary's tests in parallel — a second test here would
/// count the first one's allocations.
///
/// **The pass counts are asserted too.** Two readings that happened to run the same number of
/// passes would agree whatever the loop allocated, which is a test that cannot fail.
#[test]
fn the_loop_allocates_the_same_at_two_passes_and_at_four() {
    let profiler = dhat::Profiler::builder().testing().build();

    let mut alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
    alleles.admit(Box::from(b"T".as_slice()));
    alleles.admit(Box::from(b"C".as_slice()));

    let rows: Vec<Vec<GenericObservation>> = (0..3)
        .map(|sample| {
            (0..=sample)
                .map(|allele| observation(allele as u16, 4 + sample as u32))
                .collect()
        })
        .collect();
    let per_sample: Vec<GenericLocusSample<'_>> = rows
        .iter()
        .map(|row| GenericLocusSample {
            evidence: GenericSampleEvidence::new(row, 0.0, &[]),
            genotype_must_be_missing: false,
        })
        .collect();
    let region = GenomeRegion {
        contig: ContigId(3),
        start: Position(940),
        end: Position(940),
    };
    let evidence = LocusEvidence::generic(region, &per_sample);

    let calibration = [ReadGroupCalibration::defaulted()];
    let outbred = InbreedingF::try_new(0.0).expect("an outbred sample");
    let inbreeding = vec![outbred; 3];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    // No repeat-tract substitution rates: this fixture is a SNP/indel locus, and the map is
    // what `FrozenParameters::ssr_substitution_rate_at` answers `None` from.
    let substitution = std::collections::BTreeMap::new();
    let parameters = FrozenParameters::uncontaminated(
        &calibration,
        &inbreeding,
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
        &strata,
        &substitution,
        Ploidy::try_new(2).expect("a diploid"),
    );
    // The same run with a fraction fitted for its one library, over the default batching — one
    // batch holding all of it, which is what a run that declared nothing gets.
    let fractions = [ContaminationView {
        fraction: 0.05,
        markers_with_reads: 400,
        reads_on_markers: 1_000,
        source: ContaminationSource::ThisReadGroupsReads,
    }];
    let batching = SequencingBatches::all_together_over(1, 3);
    let contaminated = FrozenParameters::new(
        &calibration,
        &fractions,
        &batching,
        &inbreeding,
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
        &strata,
        &substitution,
        Ploidy::try_new(2).expect("a diploid"),
    );
    let arm = SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior);

    let capped_at_two = CallingLoopConfig {
        max_passes: NonZeroU32::new(2).expect("a cap of two passes"),
        ..CallingLoopConfig::DEFAULT
    }
    .validate()
    .expect("only the cap moved");
    let shipped = CallingLoopConfig::default()
        .validate()
        .expect("the shipped configuration runs");

    // **One scratch, warmed on this locus's shape before either reading.** The first call a
    // worker ever makes grows every buffer from empty, and that is an allocation of the
    // *scratch*, not of a pass — counting it would drown the number this test is about. The
    // same call warms `GenotypeTable`'s per-shape cache.
    let mut scratch = CallingScratch::<StutterSubstitutionScratch>::default();
    let _ = arm.call_locus(
        &evidence,
        &parameters,
        alleles.clone(),
        &shipped,
        &mut scratch,
    );
    // **And warmed on the contaminated shape too**, which grows three buffers the clean run
    // never sizes. Without this the contaminated arm's first reading pays for them.
    let _ = arm.call_locus(
        &evidence,
        &contaminated,
        alleles.clone(),
        &shipped,
        &mut scratch,
    );

    // **Both allele tables are cloned before either reading, and this is not a detail.** The
    // loop takes its candidates by value, and cloning one inside a measured region charges that
    // run four blocks the other never paid — which is what the first draft of this test did,
    // reporting 10 blocks against 6 and pointing the finger at the loop.
    let candidates_for_two = alleles.clone();
    let candidates_for_four = alleles.clone();
    let candidates_for_dirty_two = alleles.clone();
    let candidates_for_dirty_more = alleles;

    let before_two = dhat::HeapStats::get().total_blocks;
    let two = arm.call_locus(
        &evidence,
        &parameters,
        candidates_for_two,
        &capped_at_two,
        &mut scratch,
    );
    let after_two = dhat::HeapStats::get().total_blocks;

    let before_four = dhat::HeapStats::get().total_blocks;
    let four = arm.call_locus(
        &evidence,
        &parameters,
        candidates_for_four,
        &shipped,
        &mut scratch,
    );
    let after_four = dhat::HeapStats::get().total_blocks;

    let before_dirty_two = dhat::HeapStats::get().total_blocks;
    let dirty_two = arm.call_locus(
        &evidence,
        &contaminated,
        candidates_for_dirty_two,
        &capped_at_two,
        &mut scratch,
    );
    let after_dirty_two = dhat::HeapStats::get().total_blocks;

    let before_dirty_more = dhat::HeapStats::get().total_blocks;
    let dirty_more = arm.call_locus(
        &evidence,
        &contaminated,
        candidates_for_dirty_more,
        &shipped,
        &mut scratch,
    );
    let after_dirty_more = dhat::HeapStats::get().total_blocks;

    drop(profiler);

    assert_eq!(
        (two.passes, four.passes),
        (2, 4),
        "the two runs must differ in their pass count, or the invariant is untested"
    );
    assert_eq!(
        after_two - before_two,
        after_four - before_four,
        "two extra passes of the frequency loop allocated {} blocks against the capped run's \
         {}, so something inside a pass is allocating — the loop's buffers all belong to the \
         worker's scratch and are sized once per locus",
        after_four - before_four,
        after_two - before_two
    );

    assert!(
        dirty_more.passes > dirty_two.passes,
        "the two contaminated runs must differ in their pass count, or the invariant is \
         untested — {} against {}",
        dirty_two.passes,
        dirty_more.passes
    );
    assert_eq!(
        after_dirty_two - before_dirty_two,
        after_dirty_more - before_dirty_more,
        "{} extra passes of a contaminated locus allocated {} blocks against the capped run's \
         {}, so something inside a pass is allocating — and a contaminated pass refills three \
         buffers and assembles the whole genotype-likelihood table again",
        dirty_more.passes - dirty_two.passes,
        after_dirty_more - before_dirty_more,
        after_dirty_two - before_dirty_two
    );
}
