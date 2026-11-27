use super::*;

// ---- contamination: the mixture's second half, per locus and per sample ----

/// **A run that fitted a contamination fraction is called, and the alternative allele
/// loses ground** — which is the whole of what the correction is for: reads the neighbours
/// carry are partly explained as the neighbours' DNA rather than as this individual's
/// alternative allele (`spec/read_likelihoods.md` §3.6).
///
/// Three samples over one biallelic SNP, one batch, one library each. Two of them show
/// eight alternative reads and nothing else; the third shows twelve reference reads and two
/// alternative. **The batch is therefore mostly alternative**, so the third sample's two odd
/// reads are exactly what a contaminating neighbour would have shown — and at a 30%
/// fraction the loop attributes them there.
///
/// The same evidence is called twice, at `c = 0` and at `c = 0.3`, and what is asserted is
/// the *difference*: a single run's number says nothing about the correction, because every
/// number here also depends on the prior and on the reads.
#[test]
fn a_contaminated_run_takes_alternative_copies_off_the_cohort() {
    let (alleles, _) = generic_locus(1);
    let mostly_reference = [observation(0, 0, 12, 6, 6), observation(1, 0, 2, 1, 1)];
    let alternative_1 = [observation(1, 1, 8, 4, 4)];
    let alternative_2 = [observation(1, 2, 8, 4, 4)];
    let per_sample = [
        called_sample(&mostly_reference),
        called_sample(&alternative_1),
        called_sample(&alternative_2),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 3];
    let inbreeding = [outbred(), outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = one_batch(3, 3);

    let clean = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = worker_scratch();
    let without = shipped_arm().call_locus(
        &evidence,
        &clean,
        alleles.clone(),
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    let fractions = [contaminated_at(0.3); 3];
    let dirty = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
    let mut scratch = worker_scratch();
    let with = shipped_arm().call_locus(
        &evidence,
        &dirty,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    let clean_alternative = without.cohort_expected_copies().copies()[1];
    let dirty_alternative = with.cohort_expected_copies().copies()[1];
    assert!(
        dirty_alternative < clean_alternative,
        "attributing 30% of every library's reads to a neighbour takes alternative copies \
         off the cohort, and here it went from {clean_alternative} to {dirty_alternative}"
    );
}

/// **What a contaminated locus actually answers, pinned to the bit** — because every other
/// contaminated fixture here asserts a direction or a counter, and those pass under an
/// implementation whose per-pass assembly happens at the wrong moment.
///
/// Measured (2026-08-26): moving the head-of-pass assembly to *after* the E-step's row loop
/// — so that each pass scores against the frequencies of the pass **before** last — leaves
/// every other test in the whole `calling` module green, 733 of them, and moves this
/// fixture's cohort copies from `[3.008_738_98…, 2.991_261_01…]` to
/// `[3.008_735_69…, 2.991_264_30…]`. That is a difference of about `3.3e-6`, against a
/// tolerance here of `1e-9` — **the pass count does not move at all**, so nothing but the
/// values catches it.
///
/// The numbers are the shipped implementation's, recorded rather than hand-derived: the
/// prior goes through `lgamma`, so a hand-computed expectation would be a transcription of a
/// calculator rather than a check on this code. What makes them a check is that the two
/// implementations above disagree in the ninth digit and this notices.
#[test]
fn a_contaminated_locus_answers_the_same_numbers_it_answered_before() {
    let (alleles, _) = generic_locus(1);
    let reference = [observation(0, 0, 8, 4, 4)];
    let mixed = [observation(0, 1, 4, 2, 2), observation(1, 1, 4, 2, 2)];
    let alternative = [observation(1, 2, 8, 4, 4)];
    let per_sample = [
        called_sample(&reference),
        called_sample(&mixed),
        called_sample(&alternative),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 3];
    let fractions = [contaminated_at(0.05); 3];
    let inbreeding = [outbred(), outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = SequencingBatches::all_together_over(3, 3);
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
    let mut scratch = worker_scratch();
    let capped_alleles = alleles.clone();

    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(inference.passes, 2);
    assert!(inference.converged);
    // Re-recorded 2026-10-31, when spec §3.3's amendment began scoring a read as the average
    // over a genotype's copies, `1 − ε` from its own allele's.
    let copies = inference.cohort_expected_copies().copies();
    let expected = [3.009_985_352_332_220_6, 2.990_014_647_667_779_4];
    assert!(
        copies
            .iter()
            .zip(expected)
            .all(|(&got, want)| (got - want).abs() < 1e-9),
        "the cohort's six chromosomes at a 5% fraction: {copies:?} against {expected:?}"
    );
    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(0)]
    );
    assert_eq!(
        called(&inference, 1).0.alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    assert_eq!(
        called(&inference, 2).0.alleles(),
        [AlleleId(1), AlleleId(1)]
    );

    // **The confidences and the site quality come from the *final pass*'s table**, which is
    // the one assembled against the frequencies the loop settled on — a different table
    // from the one the last pass started from. So these are what notice the last assembly
    // going missing, where the cohort copies above cannot: they are the loop's own output
    // and it has already stopped.
    let qualities: Vec<f32> = (0..3)
        .map(|sample| called(&inference, sample).1.get())
        .collect();
    for (sample, (got, want)) in qualities
        .iter()
        .zip([20.104_8_f32, 29.945_402, 17.118_761])
        .enumerate()
    {
        assert!(
            (got - want).abs() < 1e-4,
            "sample {sample}'s confidence: {got} against {want}"
        );
    }
    // 117.716 until 2026-11-06, when the site quality stopped counting each arrangement of
    // the cohort's copies twice (`calling_quality.md` §5.2).
    assert!(
        (inference.site_quality.get() - 107.793_93_f32).abs() < 1e-3,
        "the site quality: {}",
        inference.site_quality.get()
    );

    // **And the same locus stopped at the cap, which is where the last assembly shows.**
    // The confidences and the site quality come from the *final pass*'s table, which is the
    // one assembled against the frequencies the loop settled on — a different table from the
    // one the last pass started from. At a locus that converged the two differ by less than
    // the convergence threshold, so the numbers above barely move whether or not that
    // assembly happens; stopped after one pass they have not converged, and dropping the
    // assembly moves the site quality from **119.720 to 119.742** and this sample's
    // confidence from **30.194 to 30.180** (measured 2026-08-26, before spec §3.3's
    // amendment of 2026-10-31 moved every value here; they were re-recorded then). The cohort's expected
    // copies are identical either way, which is the point: they are the loop's own output
    // and the loop has already stopped.
    let mut capped_scratch = worker_scratch();
    let capped = shipped_arm().call_locus(
        &evidence,
        &parameters,
        capped_alleles,
        &capped_at(1),
        &mut capped_scratch,
    );
    assert!(!capped.converged, "one pass is not enough for this locus");
    assert!(
        // 117.715 until 2026-11-06, as above.
        (capped.site_quality.get() - 107.792_63_f32).abs() < 1e-4,
        "the capped run's site quality: {}",
        capped.site_quality.get()
    );
    assert!(
        (called(&capped, 1).1.get() - 29.945_393_f32).abs() < 1e-4,
        "the capped run's middle sample: {}",
        called(&capped, 1).1.get()
    );
}

/// **The expensive half of the build is paid for once, whatever the pass count and
/// whether or not anything is contaminated** — and the cheap half once a pass where
/// something is.
///
/// This is the invariant `spec/calling_em_loop.md` §13's test 5 pins, in the form the
/// contamination ruling left it: `q(o)` moves with the loop, so a whole row may no longer
/// be cached, but the emissions underneath read no frequency and are still computed once
/// per `(sample, observation, candidate)` (`spec/read_likelihoods.md` §6.1).
///
/// **Two pass counts, because one cannot say "independent of the pass count".** The same
/// locus is called under a cap of 2 and then uncapped, and the assemblies are asserted as
/// `passes + 2` — one for the initialisation pass, one at the head of each pass, and one
/// against the settled frequencies before the final pass.
#[test]
fn a_contaminated_locus_builds_its_emissions_once_and_assembles_them_once_a_pass() {
    let (alleles, _) = generic_locus(2);
    // **The four-pass fixture of `the_table_is_built_once_at_a_locus_that_takes_four_passes`,
    // one library per sample** — four reads of each of three alleles at every sample, which
    // is the shape that keeps the loop moving for four passes.
    let first = [
        observation(0, 0, 4, 2, 2),
        observation(1, 0, 4, 2, 2),
        observation(2, 0, 4, 2, 2),
    ];
    let second = [
        observation(0, 1, 4, 2, 2),
        observation(1, 1, 4, 2, 2),
        observation(2, 1, 4, 2, 2),
    ];
    let third = [
        observation(0, 2, 4, 2, 2),
        observation(1, 2, 4, 2, 2),
        observation(2, 2, 4, 2, 2),
    ];
    let per_sample = [
        called_sample(&first),
        called_sample(&second),
        called_sample(&third),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 3];
    let fractions = [contaminated_at(0.05); 3];
    let inbreeding = [outbred(), outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = one_batch(3, 3);
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

    // **One scratch across both runs, which is the only shape a real worker has** — so the
    // counter's own reset is load-bearing here too.
    let mut scratch = worker_scratch();
    let mut passes_seen = Vec::new();
    for config in [capped_at(2), RunnableCallingLoopConfig::default()] {
        let inference = shipped_arm().call_locus(
            &evidence,
            &parameters,
            alleles.clone(),
            &config,
            &mut scratch,
        );
        passes_seen.push(inference.passes);
        let cost = scratch.emission_cost();
        assert_eq!(
            (
                cost.emission_builds,
                cost.emission_row_fills,
                cost.emission_evaluations
            ),
            (1, 3, 27),
            "the emissions read no frequency, so they are computed once over three rows of \
             three observations against three candidates however many passes the loop \
             takes — and it took {}",
            inference.passes
        );
        assert_eq!(
            (cost.table_assemblies, cost.row_assemblies),
            (
                u64::from(inference.passes) + 2,
                (u64::from(inference.passes) + 2) * 3
            ),
            "the table is assembled once for the initialisation pass, once at the head of \
             each of the {} passes, and once against the settled frequencies",
            inference.passes
        );
    }
    assert_eq!(
        passes_seen,
        vec![2, 7],
        "the two runs must differ in their pass count, or the invariant is untested — and \
         the uncapped one takes seven passes here against the same evidence's four with no \
         fraction fitted, because `q(o)` moves between passes as well as the frequencies"
    );
}

/// **An uncontaminated locus assembles its table once, however many passes it takes** —
/// the other half of the same claim, and the one that says the split cost a clean run
/// nothing.
///
/// The fixture is the contaminated one above with the fraction taken away, so the two can
/// be read side by side: same reads, same evidence, `table_assemblies` 1 against
/// `passes + 2`.
#[test]
fn an_uncontaminated_locus_assembles_its_table_once_however_many_passes_it_takes() {
    let (alleles, _) = generic_locus(2);
    // The same evidence as the contaminated fixture beside it, so the two can be read side
    // by side: same reads, same four passes, one fold against `passes + 2`.
    let first = [
        observation(0, 0, 4, 2, 2),
        observation(1, 0, 4, 2, 2),
        observation(2, 0, 4, 2, 2),
    ];
    let second = [
        observation(0, 1, 4, 2, 2),
        observation(1, 1, 4, 2, 2),
        observation(2, 1, 4, 2, 2),
    ];
    let third = [
        observation(0, 2, 4, 2, 2),
        observation(1, 2, 4, 2, 2),
        observation(2, 2, 4, 2, 2),
    ];
    let per_sample = [
        called_sample(&first),
        called_sample(&second),
        called_sample(&third),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 3];
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
        (
            scratch.emission_cost().table_assemblies,
            scratch.emission_cost().row_assemblies
        ),
        (1, 3),
        "with no fraction fitted the assembled row reads no frequency, so it is the same \
         value at every one of the four passes"
    );
}

/// **The initialisation assembly scores the reads alone**, because there is no estimate of
/// the contaminating population's frequencies for it to score against — no pass has run.
///
/// The assertion is that the table it produces is **the very table an uncontaminated run
/// gets**, entry for entry, which is a stronger claim than any statement about the
/// frequency buffer: it says the mixture's second term is not merely small on that
/// assembly, it is absent. A flat `q(o)`, which is what this fixture asserted before,
/// passes no part of it.
#[test]
fn the_initialisation_assembly_scores_the_reads_alone() {
    let (alleles, table) = generic_locus(2);
    let genotypes = table.view();
    let first = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let second = [observation(0, 1, 4, 2, 2), observation(1, 1, 4, 2, 2)];
    let per_sample = [called_sample(&first), called_sample(&second)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 2];
    let fractions = [contaminated_at(0.05); 2];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = one_batch(2, 2);
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

    let clean = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );

    let assembled = |parameters: &FrozenParameters<'_>,
                     frequencies: ContaminantFrequencies|
     -> Vec<LogProb> {
        let mut scratch = worker_scratch();
        scratch.prepare_for_locus(2, &alleles, &genotypes);
        scratch.claim_row_for(0, outbred());
        scratch.claim_row_for(1, outbred());
        if frequencies != ContaminantFrequencies::NoneFitted {
            scratch.prepare_contaminant_tables(parameters.batch_count(), parameters.sample_count());
        }
        fill_what_no_pass_recomputes(
            &StutterSubstitutionEmission,
            &evidence,
            parameters,
            &alleles,
            &genotypes,
            &mut scratch,
        );
        assemble_genotype_likelihood_table(
            &evidence,
            parameters,
            &genotypes,
            frequencies,
            &mut scratch,
        );
        (0..2)
            .flat_map(|row| {
                scratch
                    .sample_scoring_buffers_mut(row)
                    .genotype_likelihoods
                    .to_vec()
            })
            .collect()
    };

    let contaminated = assembled(&parameters, ContaminantFrequencies::TheReadsAlone);
    let uncontaminated = assembled(&clean, ContaminantFrequencies::NoneFitted);
    assert_eq!(
        contaminated, uncontaminated,
        "the initialisation assembly of a run whose fit found 5% contamination is the same \
         table a run that found none gets"
    );
    assert!(
        contaminated.iter().all(|value| value.0.is_finite()),
        "and the fixture actually scored something: {contaminated:?}"
    );
}

/// **A sample leaves its own copies out of its own batch**, so a contaminating read is
/// somebody else's by construction and not partly its own.
///
/// Two sequencing batches over four samples, so that the subtraction is visible against a
/// batch it did not touch. Batch 0 holds samples 0 and 1, batch 1 holds samples 2 and 3;
/// sample 0 is a reference homozygote and sample 1 an alternative one, and the two batches
/// are given different evidence so that a fill reading the wrong row is a wrong number
/// rather than the same one.
///
/// The assertion is on the **last** row's frequencies, which is what the buffer holds when
/// the fold returns — sample 3's, whose own copies come out of batch 1 and out of nothing
/// else.
#[test]
fn a_samples_own_copies_come_out_of_its_own_batch_and_no_other() {
    let (alleles, table) = generic_locus(1);
    let genotypes = table.view();
    let reference_0 = [observation(0, 0, 8, 4, 4)];
    let alternative_1 = [observation(1, 1, 8, 4, 4)];
    let reference_2 = [observation(0, 2, 8, 4, 4)];
    let reference_3 = [observation(0, 3, 8, 4, 4)];
    let per_sample = [
        called_sample(&reference_0),
        called_sample(&alternative_1),
        called_sample(&reference_2),
        called_sample(&reference_3),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 4];
    let fractions = [contaminated_at(0.05); 4];
    let inbreeding = [outbred(), outbred(), outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = two_batches_of_two();
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

    let mut scratch = worker_scratch();
    scratch.prepare_for_locus(4, &alleles, &genotypes);
    for sample in 0..4 {
        scratch.claim_row_for(sample, outbred());
    }
    scratch.prepare_contaminant_tables(parameters.batch_count(), parameters.sample_count());
    // The copies a pass would have produced, written by hand so the arithmetic below is a
    // reader's rather than the loop's.
    //
    // **Sample 3 is the alternative homozygote and no other sample's copies match its own**,
    // which is what makes the subtraction's *source* testable: with three of the four
    // samples carrying `[2, 0]`, a fill that subtracted row 0's copies from every sample
    // returned the right answer for three of them and for the asserted one, and the whole
    // suite passed under it (found by mutating `contaminant_frequency_buffers_mut(row)` to
    // `(0)`).
    let expected_copies = [[2.0, 0.0], [0.0, 2.0], [2.0, 0.0], [0.0, 2.0]];
    for (row, copies) in expected_copies.iter().enumerate() {
        let buffers = scratch.sample_scoring_buffers_mut(row);
        buffers.sample_expected_copies.copy_from_slice(copies);
    }
    fill_what_no_pass_recomputes(
        &StutterSubstitutionEmission,
        &evidence,
        &parameters,
        &alleles,
        &genotypes,
        &mut scratch,
    );
    assemble_genotype_likelihood_table(
        &evidence,
        &parameters,
        &genotypes,
        ContaminantFrequencies::TheLoopsOwnEstimate,
        &mut scratch,
    );

    assert_eq!(
        scratch.batch_allele_copies(),
        &[2.0, 2.0, 2.0, 2.0],
        "batch 0 holds one reference homozygote and one alternative homozygote, and so does \
         batch 1"
    );
    let frequencies = scratch.contaminant_allele_frequencies();
    assert_eq!(
        &frequencies[0..2],
        &[0.5, 0.5],
        "sample 3 is not in batch 0, so nothing is taken out of it"
    );
    assert_eq!(
        &frequencies[2..4],
        &[1.0, MIN_CONTAMINANT_FREQUENCY],
        "sample 3's own two alternative copies come out of batch 1, leaving sample 2's two, \
         which are all reference — the alternative is floored rather than zeroed. \
         Subtracting any other sample's copies gives a different row, which is what this \
         fixture's four distinct-enough samples are for"
    );
}

/// **A batch's copies are scattered onto the sample the row names, not onto the row's own
/// index** — and the fixture puts the uncallable sample *first*, which is where the two
/// numbers come apart.
///
/// The loop's rows are the run's sample order with the uncallable samples' gaps closed up,
/// so at a locus whose *last* sample is the uncallable one every row index equals its
/// sample's index and a scatter by row is indistinguishable from a scatter by sample.
/// Measured: with the uncallable sample last, a scatter by row passes all 4,776 tests. Here
/// sample 0 is set aside, so rows 0 and 1 are samples 1 and 2 — and a scatter by row would
/// put sample 1's copies in batch 0 and sample 2's in batch 0 as well, leaving batch 1 with
/// nothing.
#[test]
fn the_copies_are_scattered_onto_the_sample_a_row_names_and_not_onto_the_row() {
    let (alleles, table) = generic_locus(1);
    let genotypes = table.view();
    let set_aside = [observation(0, 0, 8, 4, 4)];
    let reference = [observation(0, 1, 8, 4, 4)];
    let alternative = [observation(1, 2, 8, 4, 4)];
    let per_sample = [
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(&set_aside, 0.0, &[]),
            genotype_must_be_missing: true,
        },
        called_sample(&reference),
        called_sample(&alternative),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 3];
    let fractions = [contaminated_at(0.05); 3];
    let inbreeding = [outbred(), outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    // Read groups 0 and 1 on the first plate, read group 2 on the second — so samples 0 and
    // 1 are in batch 0 and sample 2 is in batch 1.
    let groups = ReadGroups::of_libraries(&[("rg0", "s0"), ("rg1", "s1"), ("rg2", "s2")]);
    let batching = SequencingBatches::declared(
        &groups,
        &[
            std::collections::BTreeSet::from([ReadGroupId(0), ReadGroupId(1)]),
            std::collections::BTreeSet::from([ReadGroupId(2)]),
        ],
    )
    .expect("a partition of the run");
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

    let mut scratch = worker_scratch();
    scratch.prepare_for_locus(2, &alleles, &genotypes);
    scratch.claim_row_for(1, outbred());
    scratch.claim_row_for(2, outbred());
    scratch.prepare_contaminant_tables(parameters.batch_count(), parameters.sample_count());
    // Row 0 is sample 1, a reference homozygote; row 1 is sample 2, an alternative one.
    for (row, copies) in [[2.0, 0.0], [0.0, 2.0]].iter().enumerate() {
        let buffers = scratch.sample_scoring_buffers_mut(row);
        buffers.sample_expected_copies.copy_from_slice(copies);
    }
    fill_what_no_pass_recomputes(
        &StutterSubstitutionEmission,
        &evidence,
        &parameters,
        &alleles,
        &genotypes,
        &mut scratch,
    );
    assemble_genotype_likelihood_table(
        &evidence,
        &parameters,
        &genotypes,
        ContaminantFrequencies::TheLoopsOwnEstimate,
        &mut scratch,
    );

    assert_eq!(
        scratch.batch_allele_copies(),
        &[2.0, 0.0, 0.0, 2.0],
        "batch 0 holds the sample this locus set aside, which contributes nothing, and \\
         sample 1's two reference copies; batch 1 holds sample 2's two alternative ones. A \\
         scatter by row index returns [2, 2, 0, 0]"
    );
}

/// **A sample with two libraries reads the batch the *sample* ran in**, which is a different
/// axis from the batch a read group ran in even though the two are the same slice type.
///
/// At one library per sample — every sample of every benchmark cohort here — the two
/// batchings have the same length and the same entries, so handing a *sample* index to the
/// read-group batching passes every shape check. Here sample 0 has read groups 0 and 1 and
/// sample 1 has read group 2, so the sample-keyed batching is `[0, 1]` and the read-group
/// one is `[0, 0, 1]` — and reading `batch_of_each_read_group[1]` for sample 1 gives batch
/// 0, which is the batch sample 1 is not in.
///
/// **`FrozenParameters::batch_of_sample` and `batch_of_read_group` make that a type error**
/// rather than a number, since the second takes a `ReadGroupId`. This fixture is what
/// catches it if the two are ever collapsed back onto one call.
#[test]
fn a_sample_with_two_libraries_reads_its_own_samples_batch() {
    let (alleles, table) = generic_locus(1);
    let genotypes = table.view();
    let first = [observation(0, 0, 4, 2, 2), observation(1, 1, 4, 2, 2)];
    let second = [observation(0, 2, 8, 4, 4)];
    let per_sample = [called_sample(&first), called_sample(&second)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 3];
    let fractions = [contaminated_at(0.05); 3];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let groups = ReadGroups::of_libraries(&[("rg0", "s0"), ("rg1", "s0"), ("rg2", "s1")]);
    let batching = SequencingBatches::declared(
        &groups,
        &[
            std::collections::BTreeSet::from([ReadGroupId(0), ReadGroupId(1)]),
            std::collections::BTreeSet::from([ReadGroupId(2)]),
        ],
    )
    .expect("each sample's libraries ran together");
    assert_eq!(
        (
            batching.of_each_read_group().0.len(),
            batching.of_each_sample().0.len()
        ),
        (3, 2),
        "the fixture's two views are different lengths, which is what makes the mis-key \\
         reachable"
    );
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
    assert_eq!(parameters.batch_of_sample(1), BatchId(1));
    assert_eq!(parameters.batch_of_read_group(ReadGroupId(1)), BatchId(0));

    let mut scratch = worker_scratch();
    scratch.prepare_for_locus(2, &alleles, &genotypes);
    scratch.claim_row_for(0, outbred());
    scratch.claim_row_for(1, outbred());
    scratch.prepare_contaminant_tables(parameters.batch_count(), parameters.sample_count());
    for (row, copies) in [[1.0, 1.0], [0.0, 2.0]].iter().enumerate() {
        let buffers = scratch.sample_scoring_buffers_mut(row);
        buffers.sample_expected_copies.copy_from_slice(copies);
    }
    fill_what_no_pass_recomputes(
        &StutterSubstitutionEmission,
        &evidence,
        &parameters,
        &alleles,
        &genotypes,
        &mut scratch,
    );
    assemble_genotype_likelihood_table(
        &evidence,
        &parameters,
        &genotypes,
        ContaminantFrequencies::TheLoopsOwnEstimate,
        &mut scratch,
    );

    // Sample 1's row is the last written, so the buffer holds its frequencies. Its own two
    // alternative copies come out of **batch 1**, which held only them — so batch 1 falls
    // back to the reference. Taking them out of batch 0 instead would leave `[1, 0]` there
    // and a batch-1 row of `[MIN, 1]`, neither of which this asserts.
    let frequencies = scratch.contaminant_allele_frequencies();
    assert_eq!(
        &frequencies[0..2],
        &[0.5, 0.5],
        "batch 0 holds sample 0 alone, one copy of each allele, and sample 1 is not in it"
    );
    assert_eq!(
        &frequencies[2..4],
        &[1.0 - MIN_CONTAMINANT_FREQUENCY, MIN_CONTAMINANT_FREQUENCY],
        "sample 1 is alone in batch 1, so taking its own copies out leaves nothing there"
    );
}

/// **A contaminated run scored on the plain formula is refused, in release**, because it is
/// the one direction of the pairing that would be silent: the fraction the pre-pass fitted
/// goes unused at every locus of the run and nothing in the output says which formula ran.
#[test]
#[should_panic(expected = "a run that fitted no contamination is the only one")]
fn a_contaminated_run_scored_on_the_plain_formula_is_refused() {
    let (alleles, table) = generic_locus(1);
    let genotypes = table.view();
    let rows = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let fractions = [contaminated_at(0.05)];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = SequencingBatches::all_together_over(1, 1);
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

    let mut scratch = worker_scratch();
    scratch.prepare_for_locus(1, &alleles, &genotypes);
    scratch.claim_row_for(0, outbred());
    scratch.prepare_contaminant_tables(parameters.batch_count(), parameters.sample_count());
    fill_what_no_pass_recomputes(
        &StutterSubstitutionEmission,
        &evidence,
        &parameters,
        &alleles,
        &genotypes,
        &mut scratch,
    );
    assemble_genotype_likelihood_table(
        &evidence,
        &parameters,
        &genotypes,
        ContaminantFrequencies::NoneFitted,
        &mut scratch,
    );
}

/// **A cohort of one is never a contaminated run, and this is what says so end to end.**
///
/// Contamination is a comparison between samples, so a single-sample run has nothing to fit
/// a fraction from — `RunParameters::view` routes it to the uncontaminated constructor and
/// the plain formula is what runs (`spec/read_likelihoods.md` §3.6: *emit it as absent*, not
/// a fitted zero). **The single-sample case is therefore the simple case for this model, not
/// the weak one**, and the loop reaches it through the same door every locus does.
///
/// What the assertion pins is that the locus is called at all and pays for exactly one
/// assembly however many passes it takes — the same claim the three-sample fixture makes,
/// at the other end of the committed cohort range.
#[test]
fn a_cohort_of_one_is_called_on_the_plain_formula_whatever_its_fraction_would_have_been() {
    let (alleles, _) = generic_locus(1);
    let rows = [observation(0, 0, 6, 3, 3), observation(1, 0, 6, 3, 3)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());

    // No read group identified a fraction, so the run reaches the loop through the
    // uncontaminated door and carries no batching at all. (`RunParameters`' own tests are
    // where the routing from the pre-pass's maps to that door is pinned.)
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    assert!(
        parameters.contamination_is_absent() && parameters.batch_count() == 0,
        "a single-sample run carries no fraction and no batching to read one against"
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
        [AlleleId(0), AlleleId(1)],
        "six reads of each allele at one diploid sample is a heterozygote"
    );
    assert_eq!(
        (
            scratch.emission_cost().table_assemblies,
            scratch.emission_cost().row_assemblies
        ),
        (1, 1),
        "one assembly over one row, whatever the {} passes cost",
        inference.passes
    );
    assert_eq!(scratch.contaminant_batch_count(), 0);
}

/// **A sample alone in its sequencing batch has no neighbours left**, so its contaminant is
/// drawn against the reference and the floor — the conservative answer, and the right one
/// for a library nobody was sequenced beside.
///
/// Without the subtraction it would be its own contaminant: an alternative homozygote alone
/// in its batch would come back with `q(alt) = 1`, which explains its own alternative reads
/// as somebody else's.
#[test]
fn a_sample_alone_in_its_batch_is_scored_against_the_reference() {
    let (alleles, table) = generic_locus(1);
    let genotypes = table.view();
    let alternative = [observation(1, 0, 8, 4, 4)];
    let per_sample = [called_sample(&alternative), called_sample(&alternative)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [
        ReadGroupCalibration::unscaled(),
        ReadGroupCalibration::unscaled(),
    ];
    let fractions = [contaminated_at(0.05), contaminated_at(0.05)];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = one_batch_each();
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

    let mut scratch = worker_scratch();
    scratch.prepare_for_locus(2, &alleles, &genotypes);
    scratch.claim_row_for(0, outbred());
    scratch.claim_row_for(1, outbred());
    scratch.prepare_contaminant_tables(parameters.batch_count(), parameters.sample_count());
    for row in 0..2 {
        let buffers = scratch.sample_scoring_buffers_mut(row);
        buffers.sample_expected_copies.copy_from_slice(&[0.0, 2.0]);
    }
    fill_what_no_pass_recomputes(
        &StutterSubstitutionEmission,
        &evidence,
        &parameters,
        &alleles,
        &genotypes,
        &mut scratch,
    );
    assemble_genotype_likelihood_table(
        &evidence,
        &parameters,
        &genotypes,
        ContaminantFrequencies::TheLoopsOwnEstimate,
        &mut scratch,
    );

    let frequencies = scratch.contaminant_allele_frequencies();
    assert_eq!(
        &frequencies[2..4],
        &[1.0 - MIN_CONTAMINANT_FREQUENCY, MIN_CONTAMINANT_FREQUENCY],
        "sample 1 is alone in batch 1, so taking its own copies out leaves nothing and the \
         row falls back to the reference — not to `q(alt) = 1`, which is what a batch of \
         one returns without the subtraction"
    );
}

/// **A sequencing batch every one of whose samples the candidate step ruled uncallable
/// still has a row**, and it holds zeros.
///
/// The batching is the run's and the scratch rows are the locus's, so a batch can have no
/// row here at all. Summing the copies over rows would leave that batch's row of the copy
/// table unwritten, and `fill_batch_allele_copies` refuses exactly that — by name, as a
/// batching that does not describe the run. Scattering onto the run's sample axis gives it
/// zeros instead, which is what the M-step already does with an uncallable sample.
#[test]
fn a_batch_with_no_callable_sample_here_is_a_row_of_zeros_rather_than_a_refusal() {
    let (alleles, _) = generic_locus(1);
    let alternative = [observation(1, 0, 8, 4, 4)];
    let reference = [observation(0, 1, 8, 4, 4)];
    let per_sample = [
        called_sample(&alternative),
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(&reference, 0.0, &[]),
            genotype_must_be_missing: true,
        },
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [
        ReadGroupCalibration::unscaled(),
        ReadGroupCalibration::unscaled(),
    ];
    let fractions = [contaminated_at(0.05), contaminated_at(0.05)];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = one_batch_each();
    let parameters = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

    let mut scratch = worker_scratch();
    let inference = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert!(inference.per_sample[1].is_missing());
    assert_eq!(
        scratch.row_count(),
        1,
        "one row, for the one callable sample"
    );
    assert_eq!(
        &scratch.batch_allele_copies()[2..4],
        &[0.0, 0.0],
        "batch 1 holds only the sample this locus set aside, so it contributes no copies — \
         and it still has a row, which is what keeps the copy fill from refusing the run"
    );
}

/// **A worker that calls a contaminated locus and then an uncontaminated one does not carry
/// the first one's frequencies into the second.**
///
/// `prepare_for_locus` un-sizes the two contaminant tables and nothing but a contaminated
/// locus re-sizes them, so the second locus's fold reaches the uncontaminated formula
/// rather than a table of another locus's numbers. **Asserted against a fresh scratch**,
/// because "the tables are empty" is a claim about the shape and this is a claim about the
/// answer.
#[test]
fn an_uncontaminated_locus_after_a_contaminated_one_is_called_as_though_it_were_first() {
    let (alleles, _) = generic_locus(1);
    let carrier = [observation(1, 0, 8, 4, 4)];
    let reference = [observation(0, 1, 8, 4, 4)];
    let per_sample = [called_sample(&carrier), called_sample(&reference)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); 2];
    let fractions = [contaminated_at(0.3); 2];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let batching = one_batch(2, 2);
    let dirty = contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
    let clean = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );

    let mut reused = worker_scratch();
    let _ = shipped_arm().call_locus(
        &evidence,
        &dirty,
        alleles.clone(),
        &RunnableCallingLoopConfig::default(),
        &mut reused,
    );
    let after = shipped_arm().call_locus(
        &evidence,
        &clean,
        alleles.clone(),
        &RunnableCallingLoopConfig::default(),
        &mut reused,
    );
    assert_eq!(
        reused.contaminant_batch_count(),
        0,
        "the second locus fitted no fraction, so its contaminant tables were never sized"
    );

    let mut fresh = worker_scratch();
    let alone = shipped_arm().call_locus(
        &evidence,
        &clean,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut fresh,
    );
    assert_eq!(
        after.cohort_expected_copies().copies(),
        alone.cohort_expected_copies().copies(),
        "the uncontaminated locus gives the same answer on a reused worker as on a fresh one"
    );
    assert_eq!(after.passes, alone.passes);
}

/// **A sample the candidate step ruled uncallable is given no scratch row at all**, and
/// what that buys is a denominator: the cohort's expected copies are over the samples the
/// locus was called on, so they sum to the **two** chromosomes of the one called sample
/// rather than to the four of the run (`spec/calling_em_loop.md` §5.0, §9).
///
/// A version that gave the sample a row and skipped it in the E-step would leave that row
/// holding the scratch's `NaN` sentinel and the M-step would refuse the locus; one that
/// scored it would put its mass on whichever surviving genotype its reads mismatch least
/// — usually the homozygous reference — and pull the frequencies toward the reference by
/// exactly the samples carrying the rarest alleles.
#[test]
fn a_sample_the_candidate_step_ruled_uncallable_gets_no_row_and_no_vote() {
    let (alleles, _) = generic_locus(1);
    let carrier = [observation(1, 0, 8, 4, 4)];
    let loud = [observation(0, 0, 40, 20, 20)];
    let per_sample = [
        called_sample(&carrier),
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(&loud, 0.0, &[]),
            genotype_must_be_missing: true,
        },
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

    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(1), AlleleId(1)]
    );
    assert!(inference.per_sample[1].is_missing());
    assert_eq!(
        scratch.row_count(),
        1,
        "one row, for the one callable sample"
    );
    assert_eq!(
        scratch.emission_cost().emission_row_fills,
        1,
        "the uncallable sample's reads are never scored, so its emissions are never filled"
    );
    let copies: f64 = inference.cohort_expected_copies().copies().iter().sum();
    assert!(
        (copies - 2.0).abs() < 1e-12,
        "the expected copies are over the samples the locus was called on — one diploid \
         sample's two chromosomes — and they sum to {copies}"
    );
}

/// A locus at which the candidate step ruled **every** sample uncallable has nobody to
/// call, and is refused rather than emitted empty.
///
/// **The ruling that produces such a locus does not cover it.** Candidate selection cuts
/// an allele rather than refusing a locus, on the ground that most samples stay callable
/// (`candidate_alleles.md` §4.1); where none does, that argument has nothing left to rest
/// on. Nothing upstream can reach this today — the loop is not wired into a run — and the
/// refusal is loud so that the case is decided rather than discovered.
#[test]
#[should_panic(expected = "has nobody to call")]
fn a_locus_where_every_sample_was_ruled_uncallable_is_refused() {
    let (alleles, _) = generic_locus(1);
    let rows = [observation(0, 0, 4, 2, 2)];
    let set_aside = GenericLocusSample {
        evidence: GenericSampleEvidence::new(&rows, 0.0, &[]),
        genotype_must_be_missing: true,
    };
    let per_sample = [set_aside, set_aside];
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
}

/// **Four libraries over four samples, declared as two plates of two** — read groups 0 and
/// 1 on the first, 2 and 3 on the second.
///
/// The fixture that makes the two views of the batching tell each other apart is the one
/// below; this one exists so that a fill reading the wrong batch's row is a different
/// number rather than the same one.
fn two_batches_of_two() -> SequencingBatches {
    let groups =
        ReadGroups::of_libraries(&[("rg0", "s0"), ("rg1", "s1"), ("rg2", "s2"), ("rg3", "s3")]);
    SequencingBatches::declared(
        &groups,
        &[
            std::collections::BTreeSet::from([ReadGroupId(0), ReadGroupId(1)]),
            std::collections::BTreeSet::from([ReadGroupId(2), ReadGroupId(3)]),
        ],
    )
    .expect("a partition of the run")
}

/// **Two libraries over two samples, each on a plate of its own** — the batch of one, which
/// is where the leave-one-out subtraction has nothing left to leave.
fn one_batch_each() -> SequencingBatches {
    let groups = ReadGroups::of_libraries(&[("rg0", "s0"), ("rg1", "s1")]);
    SequencingBatches::declared(
        &groups,
        &[
            std::collections::BTreeSet::from([ReadGroupId(0)]),
            std::collections::BTreeSet::from([ReadGroupId(1)]),
        ],
    )
    .expect("a partition of the run")
}
