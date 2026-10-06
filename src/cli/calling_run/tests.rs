//! What the shared calling-run assembly promises, and the numbers behind the one rule here
//! that chooses rather than checks.

use super::*;

/// **What a round costs is `width × samples`, so that is what the rule holds fixed.** A
/// round holds about one observation per covered base per sample, so a width that is right
/// at sixty-three samples holds sixteen times as much at a thousand. Pinning the product
/// rather than the width is what lets one rule serve both ends of the cohort range
/// (`design_principles.md` §0).
///
/// The two clamps are what the product cannot express: below the floor the merge's own
/// default takes over, and above the ceiling there is nothing left to buy — measured on four
/// accessions over 400 kb of SL4.0, 3.29 s at 8,000 bases, 3.18 s at 32,000 and 3.22 s at
/// 64,000, the last costing 407 MB of peak resident against 340.
#[test]
fn the_round_width_holds_one_rounds_observations_to_a_budget() {
    // Between the clamps, the product is the budget.
    for samples in [40_usize, 63, 100, 500] {
        let width = u64::from(round_width_for(samples).get());
        let held = width * samples as u64;
        assert!(
            held <= u64::from(ROUND_OBSERVATION_BUDGET),
            "{samples} samples got {width} bases, holding {held} observations a round",
        );
        assert!(
            held > u64::from(ROUND_OBSERVATION_BUDGET) / 2,
            "{samples} samples got {width} bases, which leaves most of the budget unspent",
        );
    }
}

/// **A cohort big enough to be memory-bound gets the number it has today.** The floor is the
/// merge's own compiled-in default, so nothing this rule does can make a thousand-sample run
/// hold more ground than it held before the rule existed.
#[test]
fn a_large_cohort_gets_the_merges_own_default() {
    assert_eq!(
        round_width_for(2_000).get(),
        DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN,
    );
    assert_eq!(
        round_width_for(1_000).get(),
        DEFAULT_COHORT_LOCUS_BUILDER_REGIONS_LEN,
    );
}

/// **A single sample gets the ceiling, not the budget.** One sample could hold half a million
/// bases of ground within the budget; the ceiling is there because the gain has saturated
/// long before that, and because a round that wide would make the run's memory jump on the
/// one input where it is least expected.
#[test]
fn a_small_cohort_gets_the_ceiling() {
    assert_eq!(round_width_for(1).get(), WIDEST_ROUND);
    assert_eq!(round_width_for(4).get(), WIDEST_ROUND);
    // Zero files cannot reach the command, but the rule must not divide by it.
    assert_eq!(round_width_for(0).get(), WIDEST_ROUND);
}

/// **The benchmark's own cohort gets the width the sweep found.** 63 accessions over the
/// whole 8 Mb of `benchmarks/tomato1/regions.bed` ran in 193.2 s at 500 bases and 115.3 s at
/// 8,000, writing the same VCF; the rule lands on 7,936.
#[test]
fn the_tomato_cohorts_width_is_the_one_that_was_measured() {
    assert_eq!(round_width_for(63).get(), 7_936);
}

/// **A calling run writes back the SNP/indel fit's starts it read**, through the one door both
/// calling commands write their file through. A file it read carried how each start of the fit
/// behind its numbers ended; the file it writes must say the same, or every calling run would
/// strip that record from the numbers it passes on.
#[test]
fn a_calling_run_writes_back_the_fit_starts_it_read() {
    use crate::calling::parameters_file::{SnpIndelFitStart, StartOutcome};
    let read_groups = ReadGroups::of_lanes(&[("HWI.3", "TS-1", "lib3")]);
    let ploidy = Ploidy::try_new(2).expect("diploid");
    let inbreeding = DeclaredInbreeding::nothing_said();
    let starts = vec![
        SnpIndelFitStart {
            start: 1,
            ended: StartOutcome::Converged,
            passes: 41,
            agreed_with_start: None,
        },
        SnpIndelFitStart {
            start: 2,
            ended: StartOutcome::AgreedWithAnEarlierStart,
            passes: 9,
            agreed_with_start: Some(1),
        },
    ];
    let numbers = TheRunsNumbers {
        parameters: RunParameters::of_defaults(&read_groups, ploidy, &inbreeding),
        reads_behind_each_calibration: ReadsBehindEachCalibration::nothing_was_fitted(
            read_groups.len(),
        ),
        inbreeding_by_sample: inbreeding.of_each_sample(&read_groups),
        census: CensusIdentity::of_a_run_with_no_census(),
        snp_indel_fit_starts: Some(starts.clone()),
    };
    let written = numbers.parameters_file(
        &read_groups,
        &ReferenceDigest([3; 16]),
        &StrRepeatCriteria::default(),
    );
    assert_eq!(written.fitted_from.snp_indel_fit_starts, Some(starts));

    let with_none = TheRunsNumbers {
        snp_indel_fit_starts: None,
        ..numbers
    };
    let written = with_none.parameters_file(
        &read_groups,
        &ReferenceDigest([3; 16]),
        &StrRepeatCriteria::default(),
    );
    assert_eq!(written.fitted_from.snp_indel_fit_starts, None);
}
