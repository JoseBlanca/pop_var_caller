use super::*;

/// Call one repeat tract through the driver, over the two candidate lengths every tract
/// fixture here uses and against the run's fitted repeat parameters.
///
/// **The candidates and their repeat counts are supplied, not selected.** The repeat-tract
/// half of candidate selection is unwritten (`candidate_alleles_ssr.md`), so a fixture
/// states them; a later reader must not take a supplied candidate set for a step's output.
fn call_tract(
    observations_of_each_sample: &[&[SequenceObservation]],
    strata: &StratumFits,
    substitution: &std::collections::BTreeMap<
        crate::parameter_estimation::repeat_strata::StratumKey,
        crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
    >,
    calibration: &[ReadGroupCalibration],
    scratch: &mut CallingScratch<StutterSubstitutionScratch>,
) -> LocusInference {
    call_tract_over(
        observations_of_each_sample,
        &tract_alleles(),
        &tract_repeat_counts(),
        strata,
        substitution,
        calibration,
        scratch,
    )
}

/// The same, over a candidate table and repeat counts the caller chose — what the two
/// fixtures that are *about* the repeat counts need.
fn call_tract_over(
    observations_of_each_sample: &[&[SequenceObservation]],
    alleles: &CandidateAlleles,
    repeat_counts: &[NonZeroU32],
    strata: &StratumFits,
    substitution: &std::collections::BTreeMap<
        crate::parameter_estimation::repeat_strata::StratumKey,
        crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
    >,
    calibration: &[ReadGroupCalibration],
    scratch: &mut CallingScratch<StutterSubstitutionScratch>,
) -> LocusInference {
    let detail = tract_detail();
    let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
        .iter()
        .map(|observations| SsrSampleEvidence::new(observations, &detail, TEST_TRACT_LEN))
        .collect();
    let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, repeat_counts);
    let inbreeding = vec![outbred(); per_sample.len()];
    let parameters = uncontaminated_run(calibration, &inbreeding, strata, substitution, diploid());
    shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles.clone(),
        &RunnableCallingLoopConfig::default(),
        scratch,
    )
}

/// **The tolerance the existing SSR caller states**, `1e-6`
/// (`d9e7b076:src/ssr/cohort/em.rs`'s `EmCfg::tol`).
///
/// **It is that caller's number and not its rule**, and the difference is measured in
/// `a_tract_called_under_both_tolerances_gives_one_answer_and_says_what_moved`'s doc: the two
/// loops divide the same movement by different totals, so running ng at this number is
/// running ng's rule at the other caller's tolerance rather than reproducing its rule.
const SSR_CALLER_TOLERANCE: f64 = 1e-6;

/// **Three samples at four reads apiece** — one homozygous for each of the tract's two
/// lengths and one heterozygous — shared by both tolerance fixtures so that the pass counts
/// they report are comparable.
///
/// **Four reads apiece, and both neighbours were measured before settling on it.** At twelve
/// the loop settles in one pass under the shipped tolerance and two under the tighter one, so
/// there is almost nothing to compare. At two the *reads* stop deciding: the first sample is
/// called `0/1` rather than `1/1`, so the fixture would be measuring the prior rather than the
/// tolerance. Four is where the tolerances have room and the reads still settle the answer.
///
/// **Both candidates are one repeat apart and there are only two of them**, which is what a
/// tract fixture in this file is; so the two per-allele movements are equal and opposite and
/// the absolute value in the convergence test is inert here. `a_fall_larger_than_every_rise_has_not_settled`
/// is what pins that half of the rule.
fn three_samples_at_four_reads() -> (
    [SequenceObservation; 1],
    [SequenceObservation; 2],
    [SequenceObservation; 1],
) {
    (
        [tract_reads(TRACT_CANDIDATE_REPEATS[1], 4)],
        [
            tract_reads(TRACT_CANDIDATE_REPEATS[0], 2),
            tract_reads(TRACT_CANDIDATE_REPEATS[1], 2),
        ],
        [tract_reads(TRACT_CANDIDATE_REPEATS[0], 4)],
    )
}

/// **The same tract under a convergence threshold the caller chooses** — what the differential
/// against the existing caller's stopping rule needs, and nothing else uses.
///
/// Everything but the threshold is the shipped configuration, so what the two runs it is used
/// for differ in is one number.
#[allow(
    clippy::too_many_arguments,
    reason = "the same list `call_tract_over` takes, plus the one number the differential                   varies; bundling them would be a type nothing else names"
)]
fn call_tract_at_threshold(
    observations_of_each_sample: &[&[SequenceObservation]],
    strata: &StratumFits,
    substitution: &std::collections::BTreeMap<
        crate::parameter_estimation::repeat_strata::StratumKey,
        crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
    >,
    calibration: &[ReadGroupCalibration],
    convergence_threshold: f64,
    scratch: &mut CallingScratch<StutterSubstitutionScratch>,
) -> LocusInference {
    let detail = tract_detail();
    let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
        .iter()
        .map(|observations| SsrSampleEvidence::new(observations, &detail, TEST_TRACT_LEN))
        .collect();
    let repeat_counts = tract_repeat_counts();
    let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);
    let inbreeding = vec![outbred(); per_sample.len()];
    let parameters = uncontaminated_run(calibration, &inbreeding, strata, substitution, diploid());
    let config = CallingLoopConfig {
        convergence_threshold,
        ..CallingLoopConfig::DEFAULT
    }
    .validate()
    .expect("a threshold inside the configured range");
    shipped_arm().call_locus(&evidence, &parameters, tract_alleles(), &config, scratch)
}

/// **ng calls genotypes at a repeat tract**, and the three a reader can derive by hand come
/// out.
///
/// Three samples at a dinucleotide `AT` tract called over two lengths — 6 whole repeats,
/// which is the reference, and 7. The first's twenty reads all showed the 7-repeat tract;
/// the second showed ten of each; the third showed twenty of the 6-repeat one, plus two
/// reads that ran out inside the tract. At a slippage level of 4 in 100 and twenty reads
/// apiece those are not close calls, so the answers do not turn on the prior's strength:
/// **`1/1`, `0/1` and `0/0`**.
///
/// **The locus's likelihood table is built once**, whatever the pass count: a tract's row
/// reads no allele frequency where nothing is contaminated, so nothing the loop does to the
/// frequencies can move it. That is the same invariant the SNP/indel path holds, reached a
/// different way — the tract's row computes its emissions and assembles them in one call,
/// so both halves are charged together, one row of each per sample.
///
/// **The three samples hold one, two and two observations, and one of the five ran out.**
/// Both are deliberate. Equal observation counts make the emission count a three-way
/// product, in which the candidate count and the sample count can be swapped for each other
/// (`EmissionCost::emission_evaluations`' own documentation names that shape); and a fixture
/// of complete reads only makes *count the observations* and *count the complete
/// observations* the same number, where on real data over half the reads overlapping a tract
/// run out inside it.
#[test]
fn a_repeat_tract_is_called_from_its_reads() {
    let carrier = [tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)];
    let heterozygote = [
        tract_reads(TRACT_CANDIDATE_REPEATS[0], 10),
        tract_reads(TRACT_CANDIDATE_REPEATS[1], 10),
    ];
    let reference_sample = [
        tract_reads(TRACT_CANDIDATE_REPEATS[0], 20),
        tract_reads_that_ran_out(TRACT_CANDIDATE_REPEATS[0], 2),
    ];
    let mut scratch = worker_scratch();
    let inference = call_tract(
        &[&carrier, &heterozygote, &reference_sample],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut scratch,
    );

    assert_eq!(inference.region, locus_region());
    assert_eq!(inference.alleles().len(), 2, "the two supplied lengths");
    assert_eq!(
        called_alleles(&inference, 0),
        vec![1, 1],
        "twenty reads of the seven-repeat tract"
    );
    assert_eq!(
        called_alleles(&inference, 1),
        vec![0, 1],
        "ten reads of each length"
    );
    assert_eq!(
        called_alleles(&inference, 2),
        vec![0, 0],
        "twenty reads of the six-repeat tract, and two that ran out inside it"
    );
    assert!(inference.converged, "a locus this one-sided settles");
    assert_eq!(
        scratch.emission_cost(),
        EmissionCost {
            emission_builds: 1,
            emission_row_fills: 3,
            // Five observations across the three samples — 1, 2 and 2 — against two
            // candidates. Not `samples × observations × candidates`, which would be 12,
            // and not the complete observations alone, which would be 8.
            emission_evaluations: 10,
            table_assemblies: 1,
            row_assemblies: 3,
        },
        "a tract's table is built once and its two halves are charged together"
    );
}

/// **An STR bundle of one run spanning its hull is called exactly as the segment over the
/// same bases** — the same genotypes, genotype qualities, site quality, expected copies and
/// warrant — through the whole loop: candidates parsed against the run, stutter looked up
/// per run, the seed read per run, the row scored per slip set (`str_bundle_loci.md` §7).
#[test]
fn a_one_run_bundle_is_called_as_the_segment_over_its_bases() {
    use crate::locus_generation::SsrBundleDetail;
    use crate::region_typing::reference_runs::ReferenceRun;
    use crate::region_typing::run_parse::parse_into_runs;

    let carrier = [tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)];
    let heterozygote = [
        tract_reads(TRACT_CANDIDATE_REPEATS[0], 10),
        tract_reads(TRACT_CANDIDATE_REPEATS[1], 10),
    ];
    let reference_sample = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 20)];
    let observations: [&[SequenceObservation]; 3] = [&carrier, &heterozygote, &reference_sample];

    let as_segment = call_tract(
        &observations,
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut worker_scratch(),
    );

    let hull = tract_bases(TRACT_CANDIDATE_REPEATS[0]);
    let detail = SsrBundleDetail {
        left_flank: tract_detail().left_flank,
        right_flank: tract_detail().right_flank,
        runs: vec![
            ReferenceRun::new(
                0,
                hull.len() as u64,
                Motif::new(TRACT_MOTIF).expect("a dinucleotide motif"),
            )
            .expect("a run of whole copies"),
        ]
        .into_boxed_slice(),
    };
    let mut alleles = CandidateAlleles::new(
        hull.clone().into_boxed_slice(),
        LocusKind::SsrBundle(detail.clone()),
    );
    alleles.admit(tract_bases(TRACT_CANDIDATE_REPEATS[1]).into_boxed_slice());
    let parses: Vec<_> = alleles
        .iter()
        .map(|allele| parse_into_runs(allele, &hull, &detail.runs))
        .collect();
    let per_sample: Vec<SsrSampleEvidence<'_>> = observations
        .iter()
        .map(|observations| {
            SsrSampleEvidence::of_bundle(
                observations,
                crate::locus_generation::LocusLen::of_region(locus_region()),
            )
        })
        .collect();
    let evidence = LocusEvidence::ssr_bundle(locus_region(), &per_sample, &detail, &parses);
    let strata = tract_strata();
    let substitution = tract_substitution_rates();
    let calibration = tract_libraries();
    let inbreeding = vec![outbred(); per_sample.len()];
    let parameters =
        uncontaminated_run(&calibration, &inbreeding, &strata, &substitution, diploid());
    let as_bundle = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut worker_scratch(),
    );

    assert_eq!(as_bundle.per_sample, as_segment.per_sample);
    assert_eq!(
        as_bundle.uncorrected_site_quality(),
        as_segment.uncorrected_site_quality()
    );
    assert_eq!(
        as_bundle.cohort_expected_copies().copies(),
        as_segment.cohort_expected_copies().copies()
    );
    assert_eq!(as_bundle.weakest_provenance, as_segment.weakest_provenance);
    assert_eq!(as_bundle.repeat_tract, as_segment.repeat_tract);
    assert_eq!(called_alleles(&as_bundle, 1), vec![0, 1]);
}

/// **The tract's prior seed is its stratum's fitted allele length variation, restricted to the
/// lengths this locus is called over and scaled by the strength the fit holds them with** —
/// asserted as numbers, because nothing else in this file can fail if it is wrong.
///
/// The 6-repeat stratum's fitted spectrum runs `[0.10, 0.30, 0.44, 0.11, 0.05]` over offsets
/// −2 … +2 whole repeats from the reference tract length, at a concentration of 20. This
/// locus's candidates sit at offsets **0 and +1**, so they take **0.44 and 0.11**, and the
/// seed is `20 × [0.44, 0.11]` = **`[8.8, 2.2]`** — the reference length favoured four to
/// one, which is what the fit says about tracts of this stratum.
///
/// **Every genotype fixture here is decided by its reads**, deliberately, so none of them
/// moves when the seed does. Three mistakes this catches and they would not:
///
/// - the reference count taken from the wrong candidate — the last one, 7, instead of entry
///   0 — which re-centres the shape onto the 7-repeat stratum and gives a different pair
///   entirely;
/// - the spectrum read backwards, which gives `20 × [0.11, 0.44]`;
/// - the ordinary-site seed builder called instead of the tract's, which gives the run's two
///   numbers and no length in sight.
#[test]
fn the_tracts_prior_seed_is_its_fitted_allele_length_variation_over_the_candidates() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
    let mut scratch = worker_scratch();
    let _ = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(
        scratch.seed_concentration(),
        &[8.8, 2.2],
        "20 × the fitted weights at offsets 0 and +1 from the reference tract length"
    );
}

/// **A tract whose reads say nothing is decided by its fitted prior alone, and the prior it
/// is decided by is an allele length variation.**
///
/// One sample with no reads at all. Its likelihood row is flat, so its posterior is its
/// prior, and at one sample the leave-one-out subtraction removes its own vote from its own
/// prior — leaving the seed. Against the seed `[8.8, 2.2]` the three genotypes come out
/// **0.653, 0.293 and 0.053**, so the call is `0/0`.
///
/// **This is the fixture the others cannot be.** Every genotype asserted elsewhere in this
/// file rests on twelve to twenty reads a sample, where the likelihood separates the
/// genotypes by tens of nats and the largest prior shift a wrong seed can produce is about
/// four. Here the prior is the whole of it: taking the reference count from the wrong
/// candidate re-centres the shape onto the other stratum and moves every one of those three
/// numbers.
#[test]
fn a_tract_whose_reads_decide_nothing_is_decided_by_its_fitted_prior() {
    let mut scratch = worker_scratch();
    let inference = call_tract(
        &[&[]],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(
        called_alleles(&inference, 0),
        vec![0, 0],
        "the fitted spectrum favours the reference tract length three to one"
    );
}

/// **A candidate's supplied repeat count is what picks its stratum, and it is not its byte
/// length divided by the motif's.**
///
/// This is the reason `candidate_repeat_counts` travels from the locus generator through the
/// evidence and into the driver at all, and no other fixture in this file can fail if the
/// driver ignored it: everywhere else a candidate of *n* repeats spells exactly `2n` bases,
/// so counting the bases gives the right answer by accident.
///
/// **Here the two candidates spell the same number of bases and hold different numbers of
/// whole repeats.** Both are twelve bases; the first is six clean `AT` copies, the second is
/// interrupted by a `C` and holds four.
///
/// **The count is used twice and both uses are asserted**, because a driver that read it in
/// one place and derived it in the other would pass half of this:
///
/// - **the prior.** Counting bases would put both candidates at offset 0 from the reference
///   and seed them `[9.0, 9.0]`; reading the supplied counts puts the second at offset −2,
///   whose fitted weight is 0.10, and seeds them **`[9.0, 2.0]`**;
/// - **the scoring.** The run's fit has a substitution rate at the 6-repeat stratum and none
///   at the 4-repeat one, so reading the supplied counts makes the second candidate's cells
///   fall to the stated constant and the locus **`Defaulted`**. Counting bases would put
///   both candidates in the fitted stratum and claim `FittedHere`.
#[test]
fn an_interrupted_candidates_supplied_repeat_count_is_what_picks_its_stratum() {
    let mut alleles = CandidateAlleles::new(
        tract_bases(6).into_boxed_slice(),
        LocusKind::Ssr(tract_detail()),
    );
    // Twelve bases like the reference, four whole repeats because of the `C`.
    alleles.admit(Box::from(b"ATATATATCTAT".as_slice()));
    let interrupted = [
        NonZeroU32::new(6).expect("six whole repeats"),
        NonZeroU32::new(4).expect("four whole repeats"),
    ];
    assert_eq!(
        alleles
            .bases_of(AlleleId(1))
            .expect("the second candidate")
            .len(),
        alleles.reference().len(),
        "the fixture only works while the two candidates spell the same number of bases"
    );

    let reads = [tract_reads(6, 12)];
    let mut scratch = worker_scratch();
    let inference = call_tract_over(
        &[&reads],
        &alleles,
        &interrupted,
        &tract_strata(),
        // Fitted at the stratum a base count would reach, absent at the stratum the
        // supplied count reaches.
        &substitution_rates_over(&[6, 7], TRACT_READ_GROUPS),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(
        scratch.seed_concentration(),
        &[8.8, 2.0],
        "the interrupted candidate sits two whole repeats below the reference, not level \
         with it"
    );
    assert_eq!(
        inference.weakest_provenance,
        Provenance::Defaulted,
        "the interrupted candidate's stratum has no fitted substitution rate, and it is the \
         supplied count that says which stratum that is"
    );
}

/// **A library that sent this tract no read still weakens its warrant** — the one place the
/// two paths' warrant folds are not the same rule, and it is asserted rather than only
/// documented.
///
/// The run has three libraries and only the first sends a read. The fit describes that one
/// and not the other two, so the tract's parameter table — which covers every library of the
/// run, because that is the axis the row's context table is indexed on — holds four
/// defaulted cells, and the locus comes back `Defaulted`.
///
/// **It is the conservative direction and it is a departure**, from
/// `read_likelihoods.md` §4.4's *any parameter that entered the locus*: no read of libraries
/// 1 and 2 entered anything here. Narrowing it means narrowing what the parameter table
/// covers; until then this test is what stops the behaviour changing by accident.
#[test]
fn a_library_that_sent_this_tract_no_read_still_weakens_it() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
    // The same strata, but a fit that names only the library whose reads arrived.
    let describes_one_library =
        tract_strata_describing(std::collections::BTreeMap::from([(ReadGroupId(0), 0)]));

    let mut scratch = worker_scratch();
    let inference = call_tract(
        &[&reads],
        &describes_one_library,
        &substitution_rates_over(&[4, 6, 7], 1),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(
        inference.weakest_provenance,
        Provenance::Defaulted,
        "libraries 1 and 2 sent no read and the fit does not describe them, and the tract's \
         parameter table covers them anyway"
    );
}

/// **The tract's record says how many of its scoring cells fell back to a stated constant,
/// and how many of those mean the parameters and the reads came from different runs.**
///
/// **The fixture is built so that its four counts are four different numbers, and so that
/// none of them is the reference tract's repeat count**, because three are `usize`s a record
/// could carry in the wrong slot and nothing but their values would notice. Three read groups
/// and **three** candidates make **nine cells** — where two candidates would make six, which
/// is also this tract's six whole repeats, so a wrong accessor could hide behind the
/// coincidence. The run's fit then defaults a different number of them at each of the three
/// places:
///
/// - **seven cells have no fitted slippage.** The slippage fit names read group 0 only, so
///   read groups 1 and 2 contribute six; and the 6- and 7-repeat candidates' strata are
///   fitted while the 5-repeat one's is not, so read group 0 contributes one more.
/// - **six of those seven are the fit not describing this run's read groups** — the absence
///   that means the parameters and the reads came from different runs. The seventh is the
///   ordinary kind: a candidate sitting off the fitted range on perfectly good data.
/// - **three cells have no fitted substitution rate.** Those rates are fitted for read groups
///   0 and 1 at every stratum these candidates reach, so only read group 2's three fall back.
///
/// So `9, 7, 6, 3`: no two of the four coincide, none is the reference repeat count, and no
/// pair of the counts can be swapped without the assertion below changing.
#[test]
fn a_tract_reports_how_many_of_its_cells_fell_back_and_why() {
    // The reference tract at six whole repeats, one candidate a repeat shorter — whose
    // stratum this run fitted nothing at, since `tract_strata` fits 4, 6 and 7 — and one a
    // repeat longer, whose stratum it did fit. **Three candidates, not two**, so that the
    // cell count is nine rather than the six that is also this tract's repeat count.
    let mut alleles = CandidateAlleles::new(
        tract_bases(6).into_boxed_slice(),
        LocusKind::Ssr(tract_detail()),
    );
    alleles.admit(tract_bases(5).into_boxed_slice());
    alleles.admit(tract_bases(7).into_boxed_slice());
    let repeat_counts = [
        NonZeroU32::new(6).expect("six whole repeats"),
        NonZeroU32::new(5).expect("five whole repeats"),
        NonZeroU32::new(7).expect("seven whole repeats"),
    ];
    let describes_one_library =
        tract_strata_describing(std::collections::BTreeMap::from([(ReadGroupId(0), 0)]));

    let reads = [tract_reads(6, 12)];
    let mut scratch = worker_scratch();
    let inference = call_tract_over(
        &[&reads],
        &alleles,
        &repeat_counts,
        &describes_one_library,
        &substitution_rates_over(&[4, 5, 6, 7], 2),
        &tract_libraries(),
        &mut scratch,
    );

    let record = tract_record(&inference);
    assert_eq!(
        record.scoring_cells(),
        9,
        "three read groups over three candidates — and not this tract's six repeats"
    );
    assert_eq!(
        record.cells_with_no_fitted_slippage(),
        7,
        "two undescribed read groups over all three candidates, plus the 5-repeat candidate \
         under the read group the fit does describe"
    );
    assert_eq!(
        record.cells_whose_read_group_the_fit_does_not_describe(),
        6,
        "the 5-repeat candidate's own cell is the ordinary absence and is not one of these"
    );
    assert_eq!(
        record.cells_with_no_fitted_substitution_rate(),
        3,
        "read group 2's three cells; read groups 0 and 1 have a rate at every stratum here"
    );
    assert!(!record.contaminant_term_was_built());
}

/// **The two kinds of absence are counted apart, and this is the one that is ordinary.**
///
/// Every library is named by the fit, and the second candidate holds four whole repeats
/// against the reference's six — a length the run fitted a slippage stratum for and no
/// substitution rate at. So **three cells take the stated substitution rate** (one per
/// library, for that candidate) while **no cell at all is defaulted by an unknown library**,
/// where the fixture above has four of the second kind and two of the first.
///
/// Folding the two counts into one would let *the parameters and the reads came from
/// different runs* arrive looking exactly like *a candidate sits off the fitted range*.
#[test]
fn a_candidate_off_the_fit_is_not_counted_as_a_library_the_fit_never_saw() {
    let mut alleles = CandidateAlleles::new(
        tract_bases(6).into_boxed_slice(),
        LocusKind::Ssr(tract_detail()),
    );
    alleles.admit(Box::from(b"ATATATATCTAT".as_slice()));
    let interrupted = [
        NonZeroU32::new(6).expect("six whole repeats"),
        NonZeroU32::new(4).expect("four whole repeats"),
    ];

    let reads = [tract_reads(6, 12)];
    let mut scratch = worker_scratch();
    let inference = call_tract_over(
        &[&reads],
        &alleles,
        &interrupted,
        &tract_strata(),
        &substitution_rates_over(&[6, 7], TRACT_READ_GROUPS),
        &tract_libraries(),
        &mut scratch,
    );

    let record = tract_record(&inference);
    assert_eq!(record.scoring_cells(), 6);
    assert_eq!(
        record.cells_with_no_fitted_slippage(),
        0,
        "the run fitted a slippage stratum at four repeats as well as at six"
    );
    assert_eq!(record.cells_whose_read_group_the_fit_does_not_describe(), 0);
    assert_eq!(
        record.cells_with_no_fitted_substitution_rate(),
        3,
        "one per library, for the candidate whose stratum has no fitted rate"
    );
}

/// **A tract scored on a fully fitted run reports no fallback at all — including on a
/// worker that has just scored one that fell back everywhere.**
///
/// The counters live on a per-worker buffer that every tract clears and refills, so the
/// order here is the load-bearing part: **the defaulting tract runs first**. With the fully
/// fitted tract first, a counter that was never reset would read zero anyway and the test
/// would pass against the defect it exists to catch.
#[test]
fn a_second_tract_on_the_same_worker_counts_its_own_fallbacks() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
    let mut scratch = worker_scratch();

    let fell_back = call_tract(
        &[&reads],
        &tract_strata_describing(std::collections::BTreeMap::from([(ReadGroupId(0), 0)])),
        &substitution_rates_over(&[4, 6, 7], 1),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(tract_record(&fell_back).cells_with_no_fitted_slippage(), 4);

    let fitted = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut scratch,
    );
    let record = tract_record(&fitted);
    assert_eq!(record.scoring_cells(), 6);
    assert_eq!(record.cells_with_no_fitted_slippage(), 0);
    assert_eq!(record.cells_whose_read_group_the_fit_does_not_describe(), 0);
    assert_eq!(record.cells_with_no_fitted_substitution_rate(), 0);
}

/// **A SNP called on a worker that has just called a repeat tract carries no tract record**,
/// and the tract before it is what makes the assertion mean anything.
///
/// The counts come off a per-worker buffer that only a tract refills, so at an ordinary site
/// that buffer still holds the **last tract's** numbers. What decides whether a record is
/// built is therefore the locus's own prior — the tract ladder answered, or it did not —
/// and never the buffer's state. **A driver that asked the buffer instead** (*does it hold
/// any cells?*) would build a record here out of a locus already emitted, and
/// `LocusInference::new` would refuse it with a message about the wrong locus.
///
/// Every other fixture in this file gives each call a fresh worker, so on those a buffer
/// that had never held a tract answers *no cells* and the two rules agree. **This is the one
/// fixture where they can disagree**, which is why the tract's own counts are asserted
/// non-zero first: without that the worker might be clean and the test would prove nothing.
#[test]
fn a_snp_after_a_tract_on_one_worker_carries_no_tract_record() {
    let mut scratch = worker_scratch();

    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
    let tract = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(
        tract_record(&tract).scoring_cells(),
        6,
        "the worker's tract buffer is dirty, which is the whole premise of what follows"
    );

    // The same worker, now at an ordinary site.
    let (alleles, _) = generic_locus(1);
    let carrier = [observation(1, 0, 8, 4, 4)];
    let reference_sample = [observation(0, 0, 8, 4, 4)];
    let per_sample = [called_sample(&carrier), called_sample(&reference_sample)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled(); TRACT_READ_GROUPS];
    let inbreeding = [outbred(), outbred()];
    let strata = tract_strata();
    let substitution = tract_substitution_rates();
    let parameters =
        uncontaminated_run(&calibration, &inbreeding, &strata, &substitution, diploid());
    let snp = shipped_arm().call_locus(
        &evidence,
        &parameters,
        alleles,
        &RunnableCallingLoopConfig::default(),
        &mut scratch,
    );

    assert_eq!(
        snp.repeat_tract, None,
        "an ordinary site is seeded from the population curve's two moments and has no tract \
         parameters at all, whatever the worker's tract buffer still holds"
    );
}

/// **The tract's record says whether its reads were shared out with a contaminant.**
///
/// A genotype scored under a three-term mixture and one scored under two are different
/// claims about the same reads, and the call does not say which. It is a run-wide condition
/// read at the tract: the same tract, the same fit, once in a run whose parameter fit found
/// contamination and once in one that did not.
#[test]
fn a_tract_says_whether_the_contaminant_term_was_built() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 20)];

    let mut clean_scratch = worker_scratch();
    let clean = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut clean_scratch,
    );
    assert!(!tract_record(&clean).contaminant_term_was_built());

    let mut contaminated_scratch = worker_scratch();
    let contaminated = call_contaminated_tract(&[&reads], 0.08, &mut contaminated_scratch);
    assert!(tract_record(&contaminated).contaminant_term_was_built());
}

/// **The repeat-tract differential: one tract called at the tolerance the existing SSR caller
/// states and at ng's own, with the genotypes required to match and what moved reported.**
///
/// # Why a differential and not a parity oracle
///
/// The SNP/indel path has an oracle — the existing caller's own loop's recorded answers on
/// the same likelihood tables (`calling::loop_parity`). **The repeat-tract path does not**, and the
/// reason `doc/devel/ng/spec/calling_em_loop.md` §10 gives is that the two loops converge on a
/// different quantity at a different scale, so two loops stopping at different points on one
/// trajectory would disagree at any genotype near a boundary for no reason any document
/// records.
///
/// # The two stopping tests measure one movement and divide it by different totals
///
/// **This is the part §10's wording makes sound simpler than it is, and getting it wrong is
/// what an earlier version of this comment did.** Both loops take the largest per-allele
/// change in the cohort's expected allele copies between passes, and both turn it into a
/// frequency before comparing it against a tolerance. They divide by different things:
///
/// - the SSR caller adds its prior's pseudocounts to the copies before normalising
///   (`d9e7b076:src/ssr/cohort/em.rs`'s `run_pi_em`, whose `expected` **starts at**
///   `g0`), so its divisor is `chromosomes + pseudocount mass`;
/// - ng divides by the cohort's chromosomes alone.
///
/// **So at one nominal number ng's test is the stricter of the two**, by a factor of one plus
/// that mass over the chromosomes — the SSR caller sees the same movement as smaller and stops
/// sooner. Its own engine's documentation names this as a real effect rather than a rounding
/// one: a pseudocount-scaled readout *"does not feed back"*, and testing it *"let a larger
/// pseudocount damp the delta and stop the loop early"* (spec §6).
///
/// **What this fixture therefore is: ng's rule at two tolerances**, `1e-6` — the SSR caller's
/// number — and ng's own `1e-3`. **What it is not** is a reproduction of that caller's rule,
/// which would need ng's prior strength declared the counterpart of those pseudocounts; that
/// is a claim about the two models and nobody has made it. Measured, the residual is worth one
/// pass out of five here: absorbing a plausible pseudocount mass into the divisor stops the
/// tight arm at four rather than five.
///
/// # What it requires and what it reports
///
/// **Requires:** the genotypes match. That is the failing state a differential needs.
///
/// **Reports, by asserting the numbers so they cannot go stale:** the tighter tolerance takes
/// **five passes against two** on this tract — most of a run's extra work for the same
/// three genotypes — and the two stopping points land **5.1 × 10⁻⁵ of a chromosome** apart,
/// far inside the looser tolerance. *(That is how far apart the two answers
/// finished, not a promise either rule makes: a convergence rule bounds the **last step**
/// between two passes and says nothing about the distance to another rule's answer.)*
///
/// **Every number here is measured on this fixture rather than predicted from it**, which is
/// worth saying because the first draft guessed three and two at twelve reads a sample and the
/// loop settled in one.
#[test]
fn a_tract_called_under_both_tolerances_gives_one_answer_and_says_what_moved() {
    /// Measured on this fixture, not predicted from it — at the shipped outlier weight,
    /// so the counts move when that constant moves (the shipped rule took two passes at
    /// 0.05, three at 0.20, and two again at 0.07).
    const PASSES_AT_THE_TIGHT_TOLERANCE: u32 = 5;
    /// The shipped tolerance's count on the same tract.
    const PASSES_AT_THE_SHIPPED_TOLERANCE: u32 = 2;
    /// How far apart the two stopping points landed, per chromosome — see the doc above.
    const CHROMOSOMES_APART: f64 = 5.1e-5;

    let reads_of_each_sample = three_samples_at_four_reads();
    let reads: [&[SequenceObservation]; 3] = [
        &reads_of_each_sample.0,
        &reads_of_each_sample.1,
        &reads_of_each_sample.2,
    ];

    // **ng's rule at the tolerance the SSR caller states.**
    let mut tight_scratch = worker_scratch();
    let tight = call_tract_at_threshold(
        &reads,
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        SSR_CALLER_TOLERANCE,
        &mut tight_scratch,
    );

    // **The same rule at ng's own tolerance.**
    let mut shipped_scratch = worker_scratch();
    let shipped = call_tract_at_threshold(
        &reads,
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        DEFAULT_CONVERGENCE_THRESHOLD,
        &mut shipped_scratch,
    );

    // Required: the genotypes match.
    let under_each: Vec<(Vec<u16>, Vec<u16>)> = (0..3)
        .map(|sample| {
            (
                called_alleles(&tight, sample),
                called_alleles(&shipped, sample),
            )
        })
        .collect();
    for (sample, (at_tight, at_shipped)) in under_each.iter().enumerate() {
        assert_eq!(
            at_tight, at_shipped,
            "sample {sample} is called {at_tight:?} under the SSR caller's tolerance and \
             {at_shipped:?} under ng's, so this tract's call turns on the last thousandth of \
             frequency movement"
        );
    }
    assert_eq!(
        under_each
            .iter()
            .map(|(at_tight, _)| at_tight.clone())
            .collect::<Vec<_>>(),
        vec![vec![1, 1], vec![0, 1], vec![0, 0]],
        "and the answer is the one the reads give, so the fixture is comparing two settled \
         loops rather than two that both gave up"
    );
    assert!(
        tight.converged && shipped.converged,
        "a loop that hit the cap has not settled, and comparing two unsettled answers says \
         nothing about the tolerance that stopped them"
    );

    // Reported: what a thousandfold tighter tolerance costs, and how far apart it lands.
    assert_eq!(
        (tight.passes, shipped.passes),
        (
            PASSES_AT_THE_TIGHT_TOLERANCE,
            PASSES_AT_THE_SHIPPED_TOLERANCE
        ),
        "what a thousandfold tighter tolerance costs on this tract"
    );
    let moved = tight
        .cohort_expected_copies()
        .copies()
        .iter()
        .zip(shipped.cohort_expected_copies().copies())
        .map(|(at_tight, at_shipped)| (at_tight - at_shipped).abs())
        .fold(0.0_f64, f64::max);
    let chromosomes = f64::from(diploid().get()) * reads.len() as f64;
    let apart = moved / chromosomes;

    // **How far apart the two stopping points landed, asserted as a size rather than as a
    // sign.** `apart > 0` would admit a last-place float wobble as *the two runs diverged*,
    // and the whole content of this fixture is that the genotypes survive a real movement:
    // three more passes moved the frequencies somewhere, and they still called the same three
    // samples the same way.
    assert!(
        apart > CHROMOSOMES_APART / 2.0 && apart < CHROMOSOMES_APART * 2.0,
        "the two tolerances stopped {apart} of a chromosome apart, where this fixture \
         measures {CHROMOSOMES_APART}: either they now stop in the same place, and the \
         comparison is empty, or the trajectory has moved"
    );
    assert!(
        apart < DEFAULT_CONVERGENCE_THRESHOLD,
        "the two stopping points are {apart} of a chromosome apart, which is outside the \
         looser tolerance itself"
    );
}

/// **The same comparison at the tract ladder's bottom rung, where the prior's shape is flat
/// rather than fitted.**
///
/// A run whose repeat fit produced no allele length variation anywhere seeds every tract from the
/// ladder's bottom rung: a **flat** shape at one chromosome of belief. Everything else is held
/// at the fixture above's — the same three samples, the same four reads apiece, the same
/// fitted substitution rates — so the one thing that changes is the prior's shape.
///
/// **Measured: three passes against one, where the fitted spectrum takes five against two.**
/// The tighter tolerance still costs passes, and the genotypes are still the same three.
///
/// **What separates the two fixtures is the prior's shape, not its strength**, and that is
/// worth stating because the obvious guess is the other one. Swept on this tract: holding the
/// fitted shape and raising its concentration from 1 chromosome to 100 moves the pass count
/// from 5 to 4 — barely, and in the direction opposite to *more belief, more iterating*.
/// Flattening the shape at a fixed concentration drops it to 3 or 4. **An asymmetric spectrum
/// pulls the frequencies off the point the reads alone would put them, and the loop then has a
/// trajectory to walk; a flat one does not.**
///
/// *(An earlier version of this fixture changed the substitution rates as well as the shape
/// and reported one pass under both rules — the comparison empty. It was the missing rates
/// doing that, not the rung.)*
#[test]
fn the_two_tolerances_differ_at_the_ladders_bottom_rung_too() {
    /// Measured on this fixture, at ng's rule and the tolerance the SSR caller states.
    const PASSES_AT_THE_TIGHT_TOLERANCE: u32 = 3;
    /// The shipped tolerance's count on the same tract.
    const PASSES_AT_THE_SHIPPED_TOLERANCE: u32 = 1;

    let reads_of_each_sample = three_samples_at_four_reads();
    let reads: [&[SequenceObservation]; 3] = [
        &reads_of_each_sample.0,
        &reads_of_each_sample.1,
        &reads_of_each_sample.2,
    ];
    let nothing_fitted = StratumFits::over(&[], std::collections::BTreeMap::new());

    let mut tight_scratch = worker_scratch();
    let tight = call_tract_at_threshold(
        &reads,
        &nothing_fitted,
        &tract_substitution_rates(),
        &tract_libraries(),
        SSR_CALLER_TOLERANCE,
        &mut tight_scratch,
    );
    let mut shipped_scratch = worker_scratch();
    let shipped = call_tract_at_threshold(
        &reads,
        &nothing_fitted,
        &tract_substitution_rates(),
        &tract_libraries(),
        DEFAULT_CONVERGENCE_THRESHOLD,
        &mut shipped_scratch,
    );

    assert_eq!(
        rung_of(&tight),
        AlleleLengthVariationRung::StatedFlat,
        "this fixture is about the bottom rung only while it is on the bottom rung"
    );
    for sample in 0..3 {
        assert_eq!(
            called_alleles(&tight, sample),
            called_alleles(&shipped, sample),
            "sample {sample} is called differently under the two tolerances at the flattest \
             prior the ladder offers"
        );
    }
    assert!(
        tight.converged && shipped.converged,
        "a loop that hit the cap has not settled, and two unsettled answers say nothing \
         about the tolerances that stopped them"
    );
    assert_eq!(
        (tight.passes, shipped.passes),
        (
            PASSES_AT_THE_TIGHT_TOLERANCE,
            PASSES_AT_THE_SHIPPED_TOLERANCE
        ),
        "what the tighter tolerance costs where the prior's shape is flat"
    );
}

/// **A repeat tract's prior comes from its stratum's allele length variation, and the record says
/// which rung of the ladder answered.**
///
/// The run's repeat fit reached this tract's own stratum, so the rung is the stratum's own
/// fit. A run whose fit reached no stratum at all still gets a prior — the ladder always
/// answers — and the rung then says it was a stated flat shape, which is the difference
/// between a call resting on a measurement and one resting on a constant.
///
/// **A SNP/indel locus carries no rung at all**, because its prior comes from the
/// population's frequency spectrum, whose ladder has different rungs. The two are different
/// quantities and `LocusInference::new` refuses one where the other belongs.
#[test]
fn the_rung_the_tracts_prior_came_from_travels_onto_the_locus() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
    let mut scratch = worker_scratch();
    let fitted = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(rung_of(&fitted), AlleleLengthVariationRung::StratumsOwnFit);

    let mut unfitted_scratch = worker_scratch();
    let unfitted = call_tract(
        &[&reads],
        &StratumFits::over(&[], std::collections::BTreeMap::new()),
        &NO_SUBSTITUTION_RATES,
        &tract_libraries(),
        &mut unfitted_scratch,
    );
    assert_eq!(
        rung_of(&unfitted),
        AlleleLengthVariationRung::StatedFlat,
        "the ladder always answers, and says from how far down"
    );
}

/// **The same, on a contaminated run** — because the contaminant seed adds three buffers
/// that its uncontaminated twin — `a_second_tract_on_the_same_scratch_is_scored_on_its_own_parameters`,
/// further down this module — leaves empty from beginning to end, so nothing there could
/// show a clear that was dropped.
///
/// **The defaulting tract runs first here too, and it does more work in this order.** The
/// first tract's fit reaches no stratum at all, so its allele length variation is the flat one at a
/// stated concentration and its seed is spread evenly over its candidates; the second's is
/// the stratum's own fitted shape, which is not flat. A seed buffer left over from the first
/// is therefore a visibly different distribution rather than a plausible one.
#[test]
fn a_second_contaminated_tract_on_the_same_scratch_carries_none_of_the_firsts_seed() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[1], 16)];
    let defaulted_first = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 5)];

    let no_fit = StratumFits::over(&[], std::collections::BTreeMap::new());

    let mut reused = worker_scratch();
    let _ = call_contaminated_tract_under(
        &[&defaulted_first],
        0.05,
        &no_fit,
        &NO_SUBSTITUTION_RATES,
        &mut reused,
    );
    let second = call_contaminated_tract(&[&reads], 0.05, &mut reused);

    let mut fresh = worker_scratch();
    let alone = call_contaminated_tract(&[&reads], 0.05, &mut fresh);

    assert_eq!(
        second, alone,
        "the reused scratch must give the same locus, field for field and bit for bit"
    );
}

/// **The middle rung of the tract ladder reaches the record too** — a stratum the fit never
/// reached, in a run read from an older parameters file whose fit pooled that motif period's
/// tracts (no fit produces a pool now).
///
/// **It is the rung neither of the fixtures beside it can produce**, and each of them can
/// be produced by a mistake: a lookup that ignored the pool answers `StatedFlat`, and one
/// that read the pool where a stratum's own fit exists answers `PeriodsPooledTracts` at the
/// fitted tract as well. So the three rungs are asserted at three fixtures rather than
/// two (`population_diversity.md` §8, check 3).
#[test]
fn a_tract_whose_stratum_the_fit_missed_is_seeded_from_its_motif_periods_pool() {
    use crate::parameter_estimation::joint::stratum_fits::FittedAlleleLengthVariation;
    // An older parameters file whose fit reached no stratum of its own, but pooled every
    // dinucleotide tract it saw.
    let pooled = StratumFits::over(&[], std::collections::BTreeMap::new()).with_rows_of_a_file(
        [],
        std::collections::BTreeMap::from([(
            2,
            FittedAlleleLengthVariation {
                weights: vec![0.10, 0.20, 0.40, 0.20, 0.10],
                concentration: Some(12.0),
                concentration_standard_error: None,
            },
        )]),
    );

    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
    let mut scratch = worker_scratch();
    let inference = call_tract(
        &[&reads],
        &pooled,
        &NO_SUBSTITUTION_RATES,
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(
        rung_of(&inference),
        AlleleLengthVariationRung::PeriodsPooledTracts
    );
}

/// **A tract's warrant is the weakest of the repeat parameters it was scored under, and the
/// read-group calibration is not one of them.**
///
/// The two halves of that sentence are checked apart, because a fold that read the wrong
/// list would pass either one alone:
///
/// - a tract whose repeat fit reached both its strata is `FittedHere` **even though the run's
///   one library was never calibrated** — the calibration scale does not enter a tract's
///   likelihood, so charging the tract for it would report a worse warrant than the call has;
/// - a tract whose fit reached no stratum is `Defaulted` **even though the library was
///   calibrated against a measurement** — the stutter model and the substitution rate then
///   both fell to stated constants, which is what the record has to say.
#[test]
fn a_tracts_warrant_comes_from_its_repeat_parameters_and_not_from_the_calibration() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
    let calibrated = [ReadGroupCalibration {
        scale: 1.0,
        provenance: Provenance::FittedHere,
        scale_standard_error: None,
    }];

    let mut scratch = worker_scratch();
    let fitted = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut scratch,
    );
    assert_eq!(
        fitted.weakest_provenance,
        Provenance::FittedHere,
        "an uncalibrated library does not weaken a tract, whose row never reads the scale"
    );

    let mut unfitted_scratch = worker_scratch();
    let unfitted = call_tract(
        &[&reads],
        &StratumFits::over(&[], std::collections::BTreeMap::new()),
        &NO_SUBSTITUTION_RATES,
        &calibrated,
        &mut unfitted_scratch,
    );
    assert_eq!(
        unfitted.weakest_provenance,
        Provenance::Defaulted,
        "a calibrated library does not rescue a tract scored on stated constants"
    );
}

/// Call one repeat tract in a run whose fit found `fraction` of each library's reads to have
/// come from somebody else.
///
/// Every library of the run carries the same fraction, which is not what a real run looks
/// like and is not what the fixtures below are about: what they are about is what the third
/// term of the mixture does to a genotype, and one fraction is enough to show it.
fn call_contaminated_tract(
    observations_of_each_sample: &[&[SequenceObservation]],
    fraction: f64,
    scratch: &mut CallingScratch<StutterSubstitutionScratch>,
) -> LocusInference {
    call_contaminated_tract_under(
        observations_of_each_sample,
        fraction,
        &tract_strata(),
        &tract_substitution_rates(),
        scratch,
    )
}

/// The same, over a fit the caller chose — what a fixture about the seed's own buffers
/// needs, since the seed's shape is what the fit decides.
fn call_contaminated_tract_under(
    observations_of_each_sample: &[&[SequenceObservation]],
    fraction: f64,
    strata: &StratumFits,
    substitution: &std::collections::BTreeMap<
        crate::parameter_estimation::repeat_strata::StratumKey,
        crate::parameter_estimation::Estimate<crate::types::ErrorRate>,
    >,
    scratch: &mut CallingScratch<StutterSubstitutionScratch>,
) -> LocusInference {
    let detail = tract_detail();
    let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
        .iter()
        .map(|observations| SsrSampleEvidence::new(observations, &detail, TEST_TRACT_LEN))
        .collect();
    let samples = per_sample.len();
    let repeat_counts = tract_repeat_counts();
    let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);
    let calibration = tract_libraries();
    let contamination = vec![contaminated_at(fraction); TRACT_READ_GROUPS];
    let inbreeding = vec![outbred(); samples];
    let batching = one_batch(TRACT_READ_GROUPS, samples);
    let parameters = FrozenParameters::new(
        &calibration,
        &contamination,
        &batching,
        &inbreeding,
        human_like_seed(),
        strata,
        substitution,
        diploid(),
    );
    shipped_arm().call_locus(
        &evidence,
        &parameters,
        tract_alleles(),
        &RunnableCallingLoopConfig::default(),
        scratch,
    )
}

/// **A contaminated repeat tract is called, and the fitted fraction changes the answer** —
/// which is the whole of why the third term of its mixture exists.
///
/// **The failure it prevents, in one sample.** A sample carrying two copies of the 6-repeat
/// tract, in a library 8 in 100 of whose reads came from somebody else, shows twenty reads
/// of its own and four at the 7-repeat length that are not. **Without a fitted fraction the
/// model has three explanations for those four and none of them fits**: slippage to exactly
/// one repeat longer, which the stratum's fitted numbers put at about one read in two
/// hundred; the outlier term, flat over every length the tract can reach and smaller still;
/// or a second allele. It takes the third — **`0/1`**. With the fraction the pre-pass
/// measured, the contaminant explains them and the sample comes back **`0/0`**.
///
/// **The fraction's own value does the work, which is why a third run is called here.** At
/// the same four reads a fitted fraction of 1 in 100 still calls `0/1`: it is not enough
/// mass to beat a heterozygote that must also account for twenty reference reads. So this
/// fixture cannot be satisfied by a model that reads `c` as a flag.
///
/// **Where the window is, measured on this fixture at the shipped outlier weight of
/// 0.20.** Four reads is where 1 in 100 and 8 in 100 give different answers; 2 in 100
/// already recovers the homozygote. *(At the earlier weight of 0.05 the dividing line sat
/// higher — 5 in 100 still called the heterozygote — because the flat outlier floor was a
/// quarter of what it is now and the contaminant term had more work to do.)*
///
/// The other two samples are unambiguous at twenty reads and are called `0/0` and `1/1`
/// either way, so the cohort the middle sample is scored against is the same in both runs
/// and the difference is the mixture rather than the neighbours.
#[test]
fn the_fitted_fraction_stops_a_contaminants_reads_being_called_a_second_allele() {
    let reference_sample = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 20)];
    // Its own two copies, and four reads that are somebody else's.
    let contaminated_sample = [
        tract_reads(TRACT_CANDIDATE_REPEATS[0], 20),
        tract_reads(TRACT_CANDIDATE_REPEATS[1], 4),
    ];
    let carrier = [tract_reads(TRACT_CANDIDATE_REPEATS[1], 20)];
    let evidence: [&[SequenceObservation]; 3] = [&reference_sample, &contaminated_sample, &carrier];

    let mut clean_scratch = worker_scratch();
    let clean = call_tract(
        &evidence,
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut clean_scratch,
    );
    let mut contaminated_scratch = worker_scratch();
    let contaminated = call_contaminated_tract(&evidence, 0.08, &mut contaminated_scratch);
    let mut barely_scratch = worker_scratch();
    let barely = call_contaminated_tract(&evidence, 0.01, &mut barely_scratch);

    assert_eq!(
        called_alleles(&clean, 1),
        vec![0, 1],
        "with no mixture, four reads at another length are a second allele"
    );
    assert_eq!(
        called_alleles(&contaminated, 1),
        vec![0, 0],
        "with the fitted fraction they are somebody else's DNA"
    );
    // **The fraction's own value has to do work, not merely its existence.** A fixture
    // where any positive fraction gives the same answer would pass against a model that
    // read `c` as a flag; here 1 in 100 is not enough to explain four reads and 8 is.
    assert_eq!(
        called_alleles(&barely, 1),
        vec![0, 1],
        "a smaller fitted fraction cannot explain the same four reads"
    );
    for sample in [0, 2] {
        assert_eq!(
            called_alleles(&clean, sample),
            called_alleles(&contaminated, sample),
            "sample {sample} is unambiguous and is called the same either way"
        );
    }
}

/// **The contaminant seed's own numbers**, on a fitted allele length variation, over a candidate set
/// built so that no shortcut reproduces them.
///
/// **Three candidates, and each is there to kill one shortcut.** The reference spells six
/// clean `AT` repeats — twelve bases; the second spells twelve bases too but holds only four
/// whole repeats, because a `C` interrupts it; the third spells fourteen bases and seven
/// repeats.
///
/// The 6-repeat stratum's fitted spectrum is `[0.10, 0.30, 0.44, 0.11, 0.05]` over offsets
/// −2 … +2 from the reference tract length, so the three candidates sit at offsets **0, −2
/// and +1** and take raw weights **0.44, 0.10 and 0.11**, totalling 0.65. Normalised over
/// the candidate set they are `0.6769`, `0.1538` and `0.1692`; scattered onto the byte
/// lengths the reads can show, **twelve bases carries the first two summed — `0.8308` — and
/// fourteen bases carries `0.1692`**.
///
/// **What each of the three would give instead:**
///
/// - a uniform stand-in for the fit's spectrum: `2/3` and `1/3`;
/// - the shape re-centred on the last candidate rather than the reference: the 7-repeat
///   stratum's spectrum, at offsets −1, −3 and 0, one of them past the fit's reach;
/// - each candidate's share *written* at its length rather than *added* to it: twelve bases
///   keeps only `0.1538`, and the seed no longer sums to one.
///
/// **And placing a candidate by its repeat count rather than its bases cannot even
/// complete**: 6, 4 and 7 are not byte lengths this tract reaches, so the support lookup
/// finds nothing.
#[test]
fn the_contaminant_seed_is_the_fitted_spectrum_scattered_onto_the_lengths_reads_can_show() {
    let clean = tract_bases(6);
    let seven = tract_bases(7);
    // Twelve bases like the reference, four whole repeats because of the `C`.
    let interrupted = b"ATATATATCTAT".to_vec();
    let scored = [
        SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: &clean,
            repeat_count: NonZeroU32::new(6).expect("six repeats"),
        },
        SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: &interrupted,
            repeat_count: NonZeroU32::new(4).expect("four whole repeats"),
        },
        SsrCandidate {
            bundle: None,
            reference_flanks: (&[], &[]),
            flank_deletion: crate::locus_generation::FlankDeletion::NONE,
            bases: &seven,
            repeat_count: NonZeroU32::new(7).expect("seven repeats"),
        },
    ];

    let calibration = tract_libraries();
    let contamination = vec![contaminated_at(0.05); TRACT_READ_GROUPS];
    let inbreeding = [outbred()];
    let strata = tract_strata();
    let substitution = tract_substitution_rates();
    let batching = one_batch(TRACT_READ_GROUPS, 1);
    let parameters = FrozenParameters::new(
        &calibration,
        &contamination,
        &batching,
        &inbreeding,
        human_like_seed(),
        &strata,
        &substitution,
        diploid(),
    );

    let mut fits = TractScoringFits::default();
    let reference_repeats = RepeatCount(6);
    fits.gather_for_locus(
        &tract_detail().motif,
        &scored,
        TractPrior {
            reference_repeats,
            allele_length_variation: parameters.ssr_allele_length_variation_at(
                tract_detail().motif.ssr_period(),
                reference_repeats,
            ),
        },
        &parameters,
    );
    let contexts = fits.scoring_contexts(&scored);
    let locus = fits.locus_parameters(&scored, &contexts, &contamination);
    let seed = locus
        .contamination
        .expect("a run whose fit found a fraction is scored on the three-term form")
        .contaminant_length_frequencies;

    assert_eq!(seed.len(), fits.reachable_lengths().len());
    assert!(
        (seed.iter().sum::<f64>() - 1.0).abs() < 1e-12,
        "the seed is how common each length is, so it sums to one: {seed:?}"
    );
    let share_at = |length: u32| {
        let at = fits
            .reachable_lengths()
            .binary_search(&length)
            .unwrap_or_else(|_| panic!("{length} bases is not a length this tract reaches"));
        seed[at]
    };
    assert!(
        (share_at(12) - 0.54 / 0.65).abs() < 1e-12,
        "twelve bases carries both twelve-base candidates' shares summed: {}",
        share_at(12)
    );
    assert!(
        (share_at(14) - 0.11 / 0.65).abs() < 1e-12,
        "fourteen bases carries the seven-repeat candidate's alone: {}",
        share_at(14)
    );
    assert!(
        seed.iter().filter(|share| **share > 0.0).count() == 2,
        "every other reachable length carries nothing, and its reads fall to the outlier \
         floor"
    );
}

/// **A contaminated tract's genotype-likelihood table is still built once**, whatever the
/// pass count — and that is where the two paths part.
///
/// At an ordinary site the contaminant's half of the mixture is `q(o)`, the cohort's own
/// frequency for the allele an observation shows, which the loop rewrites at every pass — so
/// a contaminated SNP is assembled again at the head of each one. **At a tract the third
/// term is the fit's allele length variation for this tract's stratum, frozen before calling
/// starts**: `read_likelihoods.md` §4.5.1 weighed the cohort's own frequencies against it and
/// refused them, precisely because contamination must not move from one pass to the next. So
/// a contaminated tract costs exactly what an uncontaminated one costs.
#[test]
fn a_contaminated_tract_is_assembled_once_where_a_contaminated_snp_is_assembled_each_pass() {
    let reference_sample = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 3)];
    let heterozygote = [
        tract_reads(TRACT_CANDIDATE_REPEATS[0], 2),
        tract_reads(TRACT_CANDIDATE_REPEATS[1], 1),
    ];
    let carrier = [tract_reads(TRACT_CANDIDATE_REPEATS[1], 3)];
    let mut scratch = worker_scratch();
    let inference = call_contaminated_tract(
        &[&reference_sample, &heterozygote, &carrier],
        0.05,
        &mut scratch,
    );

    assert!(
        inference.passes > 1,
        "at three reads a sample the frequency loop has work to do, so a driver that \
         reassembled per pass would have somewhere to do it"
    );
    // **The observable is that the per-batch tables were never sized, not the assembly
    // count.** A tract's assemblies are never charged at all — the table assembly returns
    // before its own counter — so `table_assemblies` reads 1 whether the driver calls it
    // once or once a pass, and a test on that number cannot fail. What the driver genuinely
    // changed is that it treats a contaminated *tract* as a locus whose table does not move:
    // it skips the per-batch contaminant tables the reassembly would fill, and builds no
    // reassembly to fill them with.
    assert_eq!(
        scratch.contaminant_batch_count(),
        0,
        "a tract reads no per-batch contaminant frequency, so nothing sizes those tables"
    );
    assert!(
        scratch.batch_allele_copies().is_empty()
            && scratch.contaminant_allele_frequencies().is_empty(),
        "and nothing fills them"
    );
    assert_eq!(scratch.emission_cost().emission_builds, 1);
}

/// **One worker's scratch calls tract after tract, and the second is scored on its own
/// parameters** — nothing of the first survives into it.
///
/// **The defaulting tract runs first, and here is exactly what that order buys.** Run the
/// fitted tract first and every cell the second tract reads would already hold a number
/// drawn from a fit, so a stutter model or a substitution rate left over from tract one
/// still looks like a fitted one and the warrant still comes back `FittedHere`. Run the
/// *empty* fit first and every one of those cells had to be written by tract two's own
/// gather to be right at all.
///
/// **What it does not test, said so that nobody reads more into it.** Two of the buffers
/// refilled per tract are checked by an assertion rather than by this test — a stale
/// context table of the wrong width is refused by name before any row is scored — and the
/// three counters the gather resets are read by nothing outside their own module today, so
/// a dropped reset is invisible in either order. What this test genuinely covers is the
/// per-cell parameters, the tract's reachable-length support, the motif, the prior's seed
/// buffer and the row scratch's emission cache.
///
/// **The contaminant seed's three buffers are covered by
/// `a_second_contaminated_tract_on_the_same_scratch_carries_none_of_the_firsts_seed`**,
/// which is this test's shape run again on a contaminated run: this one's runs are
/// uncontaminated, so those buffers stay empty throughout it and a dropped clear could not
/// show.
///
/// It asserts the strongest form available — the second tract's whole inference, compared
/// against the same tract called on a scratch that has never seen another locus.
#[test]
fn a_second_tract_on_the_same_scratch_is_scored_on_its_own_parameters() {
    let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[1], 16)];
    let defaulted_first = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 5)];

    let mut reused = worker_scratch();
    let _ = call_tract(
        &[&defaulted_first],
        &StratumFits::over(&[], std::collections::BTreeMap::new()),
        &NO_SUBSTITUTION_RATES,
        &tract_libraries(),
        &mut reused,
    );
    let second = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut reused,
    );

    let mut fresh = worker_scratch();
    let alone = call_tract(
        &[&reads],
        &tract_strata(),
        &tract_substitution_rates(),
        &tract_libraries(),
        &mut fresh,
    );

    assert_eq!(
        second, alone,
        "the reused scratch must give the same locus, field for field and bit for bit"
    );
    assert_eq!(
        reused.emission_cost(),
        fresh.emission_cost(),
        "the cost counters are reset per locus, so the second tract's are the first's"
    );
}
