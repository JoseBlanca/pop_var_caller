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
        fitted_density: None,
        duplicated_positions: None,
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

/// **A calling run writes back both informational tables it read** (format 4,
/// `parameter_prepass_joint_fit.md` §2.1.3, traps): what the density determined and the duplicated
/// class's numbers, through the same door as the fit's starts. No run reads them, so a writer that
/// forgot them would compile and every calling run would strip them from the file it passes on.
#[test]
fn a_calling_run_writes_back_the_informational_tables_it_read() {
    use crate::calling::parameters_file::{
        FittedDensity, FittedDuplicatedPositions, InformationalNumber, PairedNumber,
    };
    let read_groups = ReadGroups::of_lanes(&[("HWI.3", "TS-1", "lib3")]);
    let ploidy = Ploidy::try_new(2).expect("diploid");
    let inbreeding = DeclaredInbreeding::nothing_said();
    let number = |value: f64| InformationalNumber {
        value,
        standard_error: Some(value / 10.0),
    };
    let marked = |value: f64| PairedNumber {
        value,
        standard_error: None,
        separately_determined: false,
    };
    let density = FittedDensity {
        stated_frequency: 0.01,
        invariant_share: marked(0.84),
        fixed_non_reference_share: number(3e-5),
        segregating_shape_a: marked(0.02),
        segregating_shape_b: number(4.28),
        share_segregating_above_stated_frequency: number(8.4e-3),
        expected_heterozygosity: number(1.2e-3),
        expected_alternative_frequency: number(7.7e-4),
    };
    let duplicated = FittedDuplicatedPositions {
        stated_frequency: 0.01,
        duplicated_share: marked(0.015),
        carrier_shape_a: marked(0.02),
        carrier_shape_b: number(1.41),
        share_with_carriers_above_stated_frequency: number(1.2e-3),
    };
    let numbers = TheRunsNumbers {
        parameters: RunParameters::of_defaults(&read_groups, ploidy, &inbreeding),
        reads_behind_each_calibration: ReadsBehindEachCalibration::nothing_was_fitted(
            read_groups.len(),
        ),
        inbreeding_by_sample: inbreeding.of_each_sample(&read_groups),
        census: CensusIdentity::of_a_run_with_no_census(),
        snp_indel_fit_starts: None,
        fitted_density: Some(density.clone()),
        duplicated_positions: Some(duplicated.clone()),
    };
    let written = numbers.parameters_file(
        &read_groups,
        &ReferenceDigest([3; 16]),
        &StrRepeatCriteria::default(),
    );
    assert_eq!(written.ordinary_site_prior.fitted_density, Some(density));
    assert_eq!(written.duplicated_positions, Some(duplicated));
}

/// **`--min-site-quality` takes zero and any finite positive quality, and refuses the rest** —
/// zero is how a user turns the threshold off, so it must parse.
#[test]
fn the_min_site_quality_parser_takes_zero_and_refuses_a_negative_or_non_number() {
    assert_eq!(parse_min_site_quality("0"), Ok(0.0));
    assert_eq!(parse_min_site_quality("1"), Ok(1.0));
    assert_eq!(parse_min_site_quality("30.5"), Ok(30.5));
    for refused in ["-1", "inf", "NaN", "high"] {
        assert!(
            parse_min_site_quality(refused).is_err(),
            "`{refused}` must be refused"
        );
    }
}

/// **The scratch directory is the run's own**: made inside `--tmp-dir` when one is given and
/// inside the VCF's directory when not, named `pvc_tmp_…`, and gone once the run lets it go —
/// while the directory it was made in stays.
#[test]
fn the_scratch_directory_is_made_where_asked_and_removed_after() {
    std::fs::create_dir_all("tmp").expect("the project's tmp directory");
    let parent = tempfile::tempdir_in("tmp").expect("a directory to work in");
    let output = parent.path().join("out").join("calls.vcf.gz");
    std::fs::create_dir_all(output.parent().expect("a parent")).expect("the output directory");

    let beside = scratch_directory_for(&output, None).expect("made beside the VCF");
    assert_eq!(beside.path().parent(), output.parent());
    let given = parent.path().join("scratch");
    std::fs::create_dir_all(&given).expect("the --tmp-dir");
    let inside = scratch_directory_for(&output, Some(&given)).expect("made inside --tmp-dir");
    assert_eq!(inside.path().parent(), Some(given.as_path()));
    assert!(
        inside
            .path()
            .file_name()
            .is_some_and(|name| name.to_string_lossy().starts_with("pvc_tmp_"))
    );

    let (made_beside, made_inside) = (beside.path().to_path_buf(), inside.path().to_path_buf());
    drop((beside, inside));
    assert!(!made_beside.exists() && !made_inside.exists());
    assert!(given.is_dir(), "the user's --tmp-dir itself is left alone");
}

/// A `--tmp-dir` that is not there is refused, naming it.
#[test]
fn a_tmp_dir_that_is_not_there_is_refused() {
    let refused = scratch_directory_for(Path::new("calls.vcf"), Some(Path::new("no/such/dir")))
        .expect_err("there is nowhere to make it");
    assert!(matches!(refused, CallingRunError::ScratchDirectory { .. }));
}

#[test]
fn standard_output_is_named_as_such() {
    assert_eq!(where_the_calls_went(Path::new("-")), "standard output");
    assert_eq!(where_the_calls_went(Path::new("calls.vcf")), "calls.vcf");
}
