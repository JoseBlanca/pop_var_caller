use super::*;

// ────────────────────────────────────────────────────────────────────────────────
// The driver: the table built once, the outer rounds structurally off
// ────────────────────────────────────────────────────────────────────────────────

/// **The driver calls genotypes from evidence** — reads in, `LocusInference` out, with the
/// likelihood table built from the reads rather than handed in.
///
/// Two samples at a biallelic SNP: the first shows 8 reads of the alternative and none of
/// the reference, the second 8 of the reference and none of the alternative. The reads are
/// one-sided enough that the calls do not depend on the prior's strength — `1/1` and
/// `0/0` — and the cohort's expected copies come out **`[2.0151, 1.9849]`**, which is the
/// two calls' four chromosomes split very nearly evenly, as the two samples' reads split
/// them. The 0.015 of a copy away from an even split is the seed's pull toward the
/// reference, which at eight reads a sample is most of what the prior is still worth.
///
/// **The table is built once**: `EmissionCost` reports one table build over two row
/// builds, whatever the pass count, which is the property the driver exists to keep.
#[test]
fn the_driver_calls_genotypes_from_reads_and_builds_the_table_once() {
    let (alleles, _) = generic_locus(1);
    let carrier = [observation(1, 0, 8, 4, 4)];
    let reference_sample = [observation(0, 0, 8, 4, 4)];
    let per_sample = [called_sample(&carrier), called_sample(&reference_sample)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(1), AlleleId(1)]
    );
    assert_eq!(
        called(&inference, 1).0.alleles(),
        [AlleleId(0), AlleleId(0)]
    );
    assert!(inference.converged);
    assert_eq!(inference.passes, 2);
    let copies = inference.cohort_expected_copies().copies();
    assert!(
        (copies[0] - 2.017_218_810_602_41).abs() < 1e-12
            && (copies[1] - 1.982_781_189_397_590_2).abs() < 1e-12,
        "the cohort's four chromosomes split as the reads split them: {copies:?}"
    );
    assert_eq!(
        scratch.emission_cost(),
        EmissionCost {
            emission_builds: 1,
            emission_row_fills: 2,
            emission_evaluations: 4,
            table_assemblies: 1,
            row_assemblies: 2,
        },
        "one emission build over two samples, each of one observation against two \
         candidates, and one fold of the table — nothing here is contaminated, so the \
         folded row reads no frequency and does not move. `passes` was {}",
        inference.passes
    );
}

/// **The arm says which arm it is**, because the seam exists to compare three of them and
/// a result that cannot name its own is not auditable.
#[test]
fn the_shipped_arm_names_itself_and_is_object_safe() {
    let arm: Box<dyn LocusGenotyper<StutterSubstitutionScratch>> = Box::new(shipped_arm());
    assert!(arm.name().contains("summarise"));
}

/// **A locus whose alternatives are the reference's own length is a substitution, and one
/// carrying a different length is an insertion or a deletion.**
///
/// The two classes take the same seed today, so this moves no number — which is exactly
/// why it needs a test: nothing downstream would notice it being wrong until the split
/// arrives (`spec/calling_priors.md` §4.2).
#[test]
fn the_variant_class_is_read_off_the_candidate_lengths() {
    let (substitution, _) = generic_locus(1);
    assert_eq!(variant_class_of(&substitution), VariantClass::Substitution);

    let mut deletion = CandidateAlleles::new(Box::from(b"ACGT".as_slice()), LocusKind::Generic);
    deletion.admit(Box::from(b"A".as_slice()));
    assert_eq!(
        variant_class_of(&deletion),
        VariantClass::InsertionOrDeletion
    );

    let reference_only = CandidateAlleles::new(Box::from(b"ACGT".as_slice()), LocusKind::Generic);
    assert_eq!(
        variant_class_of(&reference_only),
        VariantClass::Substitution,
        "a locus called over the reference alone carries no length change"
    );
}

/// **The gap comes first, which is the shape that separates a row from its sample.**
///
/// Every other set-aside fixture in this file puts the uncallable sample **last**, where a
/// row's index and its sample's index are the same number and a table filled by row gives
/// the right answer by accident. Here sample 0 is set aside holding 8 reference reads and
/// sample 1 is called holding 8 alternative ones, so the one scratch row is row 0 and it
/// belongs to run sample 1.
///
/// **Measured:** a build that read `per_sample[row]` calls the surviving sample `0/0` on
/// reads it never saw, with copies `[2.0000, 3.9e-6]` against this fixture's
/// `[0.0077, 1.9923]` — a systematic permutation, and no panic anywhere.
#[test]
fn call_locus_scores_each_row_against_the_sample_that_claimed_it_when_the_gap_comes_first() {
    let (alleles, _) = generic_locus(1);
    let reference_reads = [observation(0, 0, 8, 4, 4)];
    let alternative_reads = [observation(1, 0, 8, 4, 4)];
    let per_sample = [
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(&reference_reads, 0.0, &[]),
            genotype_must_be_missing: true,
        },
        called_sample(&alternative_reads),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert!(inference.per_sample[0].is_missing());
    assert_eq!(
        called(&inference, 1).0.alleles(),
        [AlleleId(1), AlleleId(1)],
        "row 0 belongs to run sample 1, whose reads are all alternative"
    );
    let copies = inference.cohort_expected_copies().copies();
    assert!(
        (copies[0] - 0.008_884_578_139_146_705).abs() < 1e-12
            && (copies[1] - 1.991_115_421_860_853_3).abs() < 1e-12,
        "the one called sample's two chromosomes are nearly all alternative: {copies:?}"
    );
}

/// **A record may not claim its parameters were fitted when they were defaulted.**
///
/// The fixture's calibration is `ReadGroupCalibration::unscaled()`, whose own provenance
/// is `Defaulted` — nothing could be fitted and a stated constant was used — and the reads
/// scored at this locus all came from that read group. So the weakest warrant behind the
/// call is `Defaulted`, and a consumer that treated it as fitted is the failure
/// `LocusInference::weakest_provenance` exists to prevent.
#[test]
fn call_locus_reports_the_weakest_warrant_of_the_parameters_that_reached_the_locus() {
    let (alleles, _) = generic_locus(1);
    let het = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let per_sample = [called_sample(&het)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    assert_eq!(calibration[0].provenance, Provenance::Defaulted);
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(inference.weakest_provenance, Provenance::Defaulted);
}

/// **A locus no read reached claims the strongest warrant**, because nothing weaker
/// entered it: every sample is decided by the prior alone, and the prior's own seed
/// carries no provenance to report. The fold's identity, not a claim about a measurement.
#[test]
fn a_locus_no_read_reached_has_no_weaker_warrant_to_report() {
    let (alleles, _) = generic_locus(1);
    let per_sample = [GenericLocusSample {
        evidence: GenericSampleEvidence::empty(),
        genotype_must_be_missing: false,
    }];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(inference.weakest_provenance, Provenance::FittedHere);
    assert!(inference.artifact_test_counts().is_none());
}

/// **The run's configuration reaches the loop**, and a cap is the one setting that shows
/// it: the same fixture settles on its second pass at the default cap, so a cap of one
/// stops it short and the record has to say so.
///
/// **Measured:** with the config dropped and the default used instead, this locus comes
/// back `passes = 2, converged = true` — a stronger claim than the loop earned.
#[test]
fn call_locus_honours_the_runs_pass_cap() {
    let (alleles, _) = generic_locus(1);
    let carrier = [observation(1, 0, 8, 4, 4)];
    let reference_sample = [observation(0, 0, 8, 4, 4)];
    let per_sample = [called_sample(&carrier), called_sample(&reference_sample)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference =
        shipped_arm().call_locus(&evidence, &parameters, alleles, &capped_at(1), &mut scratch);

    assert_eq!(inference.passes, 1);
    assert!(
        !inference.converged,
        "the cap stopped this locus, and the record has to say so"
    );
}

/// **Two coefficients that differ, because every other driver fixture is outbred.**
///
/// The loop reads a row's inbreeding coefficient off the scratch and the final pass reads
/// the same sample's off the run, so a compaction indexed wrongly makes the two disagree
/// about one sample with neither of them noticing. Both samples here show the same
/// balanced reads, so the only thing separating their genotype qualities is the
/// coefficient each was scored against.
#[test]
fn call_locus_claims_each_row_with_its_own_samples_inbreeding_coefficient() {
    let (alleles, _) = generic_locus(1);
    let balanced = [observation(0, 0, 3, 2, 1), observation(1, 0, 3, 1, 2)];
    let per_sample = [called_sample(&balanced), called_sample(&balanced)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [
        InbreedingF::try_new(0.0).expect("an outbred sample"),
        InbreedingF::try_new(0.9).expect("a highly inbred sample"),
    ];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    let outbred_quality = called(&inference, 0).1.get();
    let inbred_quality = called(&inference, 1).1.get();
    assert!(
        (outbred_quality - 32.109_31).abs() < 1e-3 && (inbred_quality - 20.539_87).abs() < 1e-3,
        "the same reads under F = 0 and F = 0.9: {outbred_quality} and {inbred_quality}"
    );
    let copies = inference.cohort_expected_copies().copies();
    assert!(
        (copies[0] - 2.003_347_292_611_038_3).abs() < 1e-12,
        "scoring both rows at F = 0 gives a different split: {copies:?}"
    );
}

/// **Unequal observation counts, which is the shape spec §13's test 5 asks for.**
///
/// A sample of one observation beside a sample of two: `Σ_s obs_s × candidates` is
/// `(1 + 2) × 2 = 6`, where a version charging the first row's count for every row would
/// report `2 × 1 × 2 = 4`. A fixture whose samples hold equal counts cannot tell the two
/// apart, which is what the counter's own documentation says and what the driver's first
/// fixture does.
#[test]
fn emission_evaluations_sum_each_samples_own_observation_count() {
    let (alleles, _) = generic_locus(1);
    let thin = [observation(1, 0, 8, 4, 4)];
    let thick = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let per_sample = [called_sample(&thin), called_sample(&thick)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let _ = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(
        scratch.emission_cost(),
        EmissionCost {
            emission_builds: 1,
            emission_row_fills: 2,
            emission_evaluations: 6,
            table_assemblies: 1,
            row_assemblies: 2,
        }
    );
}

/// **A read that saw only part of the locus is still an emission the builder was asked
/// for**, so it is charged: two whole-span observations and one partial, over two
/// candidates, is 6 — dropping the partials from the charge reports 4.
#[test]
fn emission_evaluations_charge_the_partial_observations_too() {
    let (alleles, _) = generic_locus(1);
    let rows = [observation(0, 0, 2, 1, 1), observation(1, 0, 6, 5, 2)];
    let partials = [PartialObservation {
        held_flank: None,
        witnessed_in_locus: WitnessedLocusPositions::from_half_open_runs([(0_u16, 1_u16)])
            .expect("one witnessed position"),
        read_group: ReadGroupId(0),
        bases: Box::from(b"A".as_slice()),
        num_reads: 50,
        q_sum: -100.0,
        sequence_may_run_on_past: None,
    }];
    let per_sample = [GenericLocusSample {
        evidence: GenericSampleEvidence::new(&rows, 0.0, &partials),
        genotype_must_be_missing: false,
    }];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let _ = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(
        scratch.emission_cost(),
        EmissionCost {
            emission_builds: 1,
            emission_row_fills: 1,
            emission_evaluations: 6,
            table_assemblies: 1,
            row_assemblies: 1,
        },
        "two whole-span observations and one partial, over two candidates"
    );
}

/// **The same one build at four passes as at two** — which is what "independent of the
/// pass count" means, and what a single fixture at two passes cannot say. Three samples
/// over three alleles, each showing one read of each: the loop takes four passes and the
/// table is still built once.
///
/// This is also the driver's only fixture above two samples and above two alleles.
#[test]
fn the_table_is_built_once_at_a_locus_that_takes_four_passes() {
    let (alleles, _) = generic_locus(2);
    let mixed = [
        observation(0, 0, 4, 2, 2),
        observation(1, 0, 4, 2, 2),
        observation(2, 0, 4, 2, 2),
    ];
    let per_sample = [
        called_sample(&mixed),
        called_sample(&mixed),
        called_sample(&mixed),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred(), outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(inference.passes, 4);
    assert_eq!(
        scratch.emission_cost(),
        EmissionCost {
            emission_builds: 1,
            emission_row_fills: 3,
            emission_evaluations: 27,
            table_assemblies: 1,
            row_assemblies: 3,
        },
        "four passes, one emission build and one fold — three rows of three observations \
         against three candidates"
    );
    let total: f64 = inference.cohort_expected_copies().copies().iter().sum();
    assert!(
        (total - 6.0).abs() < 1e-12,
        "three diploid samples carry six chromosomes, and they are all accounted for"
    );
}

/// **A cohort of one, called end to end** — the hardest corner of the range this caller
/// commits to, and the only cohort size at which the prior's leave-one-out subtraction is
/// between a number and itself. The concentration comes back as the seed, so the loop
/// reaches its fixed point by arithmetic on the second pass.
#[test]
fn call_locus_calls_a_cohort_of_one() {
    let (alleles, _) = generic_locus(1);
    let het = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let per_sample = [called_sample(&het)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    assert_eq!(inference.passes, 2);
    assert!(inference.converged);
    let copies = inference.cohort_expected_copies().copies();
    assert!(
        (copies[0] - 1.020_306_935_702_036).abs() < 1e-12
            && (copies[1] - 0.979_693_064_297_963_8).abs() < 1e-12,
        "one diploid sample's two chromosomes, split as its reads split them: {copies:?}"
    );
}

/// **Each path's evidence accessor refuses the other path's locus**, rather than reading a
/// field that happens to be there.
///
/// Neither is reachable through the driver, which branches on the evidence before reading
/// it — that is what makes them worth asserting: a branch dropped later surfaces here,
/// naming the path, rather than as a genotype scored under the wrong read model.
#[test]
#[should_panic(expected = "reached the SNP/indel path's evidence")]
fn a_tract_handed_to_the_snp_paths_evidence_accessor_is_refused() {
    let detail = tract_detail();
    let per_sample = [SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN)];
    let repeat_counts = tract_repeat_counts();
    let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);
    let _ = generic_evidence_of(&evidence);
}

/// The mirror of it, on the other path.
#[test]
#[should_panic(expected = "reached the repeat-tract path's evidence")]
fn a_snp_handed_to_the_tract_paths_evidence_accessor_is_refused() {
    let rows = [observation(0, 0, 4, 2, 2)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let _ = tract_evidence_of(&evidence);
}

/// **Spec §5.0.1's ruling, as a unit**: a discovery round at a tract can put back a length
/// the cap cut, so no sample is locked out of the locus for the rest of its calling.
///
/// **Asked here rather than only through a called tract**, because a tract that came back
/// with every sample called is consistent both with this rule and with a tract where every
/// sample happened to be callable anyway — the two are the same output.
#[test]
fn is_callable_rules_no_sample_out_on_a_repeat_tract() {
    let detail = SsrDetail {
        motif: Motif::new(b"AT").expect("a dinucleotide motif"),
        left_flank: Box::from(b"CCCGGG".as_slice()),
        right_flank: Box::from(b"TTTAAA".as_slice()),
    };
    let per_sample = [
        SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN),
        SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN),
    ];
    let repeat_counts = tract_repeat_counts();
    let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);

    assert!(is_callable(&evidence, 0));
    assert!(is_callable(&evidence, 1));
}
