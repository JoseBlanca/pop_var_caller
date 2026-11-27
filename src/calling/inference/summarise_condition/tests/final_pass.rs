use super::*;

// ────────────────────────────────────────────────────────────────────────────────
// C3b — the final pass: the calls, the site quality, and the artifact summary
// ────────────────────────────────────────────────────────────────────────────────

/// **One sample's call, its confidence, and its artifact counts, on numbers a reader can
/// follow** — the final pass's own hand-computed case.
///
/// The likelihoods and log-priors are the ones
/// [`one_samples_e_step_matches_the_arithmetic_done_by_hand`] uses, so the posterior is
/// the same `1/7, 2/7, 4/7`: the winner is genotype 2, which at a diploid biallelic
/// locus is `1/1`, and it takes `4/7` of the probability.
///
/// - **the genotype quality** is `−10·log₁₀(1 − 4/7)` = `10·log₁₀(7/3)` = **3.6798**,
///   which is what comes back to the last digit an `f32` holds: `3.6797678`;
/// - **the call** is two copies of allele 1;
/// - **the artifact counts** are the sample's own rows: 3 reference reads of which 2 are
///   forward and 1 placed left, 6 alternative reads of which 5 are forward and 2 placed
///   left, 9 in total — and because the call carries both copies of the alternative, the
///   genotypes expect all 9 of those reads to carry it.
///
/// **The reference row's forward and placed-left counts are deliberately different.** They
/// feed two different binomial tests — strand bias and read-position bias — and a fixture
/// in which they are equal cannot tell the two apart: swapping the two `+=` lines that
/// accumulate them leaves such a suite green.
///
/// **A fixed prior rather than a shipped one**, so every intermediate stays hand-checkable:
/// both shipped priors go through `lgamma`, and a hand-computed expectation for one of
/// them would be a transcription of a calculator rather than a check on this pass.
#[test]
fn one_samples_call_and_its_confidence_match_the_arithmetic_done_by_hand() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let two = std::f64::consts::LN_2;
    let rows = [observation(0, 0, 3, 2, 1), observation(1, 0, 6, 5, 2)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![0.0, two, 0.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &FixedLogPriors(vec![0.0, 0.0, 2.0 * two]),
    );

    let (genotype, quality) = called(&inference, 0);
    assert_eq!(genotype.alleles(), [AlleleId(1), AlleleId(1)]);
    let by_hand = 10.0 * float::log10(7.0 / 3.0);
    assert!(
        (f64::from(quality.get()) - by_hand).abs() < 1e-4,
        "the genotype quality is {} and the arithmetic gives {by_hand}",
        quality.get()
    );

    let counts = inference
        .artifact_test_counts()
        .expect("six reads reached the alternative allele");
    assert_eq!(counts.primary_alternative, AlleleId(1));
    assert_eq!(counts.reference_reads, 3.0);
    assert_eq!(counts.reference_forward_reads, 2.0);
    assert_eq!(counts.reference_placed_left_reads, 1.0);
    assert_eq!(counts.alternative_reads, 6.0);
    assert_eq!(counts.alternative_forward_reads, 5.0);
    assert_eq!(counts.alternative_placed_left_reads, 2.0);
    assert_eq!(counts.total_reads, 9.0);
    assert_eq!(counts.genotype_expected_alternative_reads, 9.0);
}

/// **The primary alternative is the allele the most reads reached, and the fixture makes
/// it the second one.**
///
/// Every other fixture in this module is biallelic, where allele 1 is the only
/// alternative there is and a hard-coded `AlleleId(1)` would pass every one of them. Here
/// allele 2 draws 9 reads across the cohort against allele 1's 6, and the first sample's
/// counts are pooled from **two read groups**, which a fixture with one row per allele
/// would also let a reader-of-the-first-row pass.
///
/// Both samples are called homozygous for allele 2, so preferring a called allele and
/// preferring the most reads agree here; where they disagree is
/// `the_primary_alternative_is_one_a_call_carries_before_one_with_more_reads`.
#[test]
fn the_primary_alternative_is_the_allele_the_most_reads_reached_over_two_read_groups() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    assert_eq!(view.allele_count(), 3);
    let first = [
        observation(0, 0, 1, 1, 0),
        observation(0, 1, 1, 0, 1),
        observation(1, 0, 2, 1, 1),
        observation(1, 1, 2, 1, 1),
        observation(2, 0, 3, 2, 1),
        observation(2, 1, 2, 1, 0),
    ];
    let second = [observation(1, 0, 2, 1, 1), observation(2, 0, 4, 3, 2)];
    let per_sample = [called_sample(&first), called_sample(&second)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    // Six genotypes at three alleles; both samples' reads favour the one carrying
    // allele 2 twice, which is the last of them.
    let favouring_two_twos = vec![-6.0, -6.0, -6.0, -6.0, -6.0, 0.0];
    let inference = call_locus(
        &[favouring_two_twos.clone(), favouring_two_twos],
        &evidence,
        &[1.0, 0.5, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    let counts = inference
        .artifact_test_counts()
        .expect("both alternatives drew reads");
    assert_eq!(
        counts.primary_alternative,
        AlleleId(2),
        "allele 2 drew 9 reads across the cohort and allele 1 drew 6"
    );
    // Allele 2's rows only: 3 + 2 in the first sample, 4 in the second.
    assert_eq!(counts.alternative_reads, 9.0);
    assert_eq!(counts.alternative_forward_reads, 6.0);
    assert_eq!(counts.alternative_placed_left_reads, 3.0);
    // The reference's two rows, which is the pooling across read groups.
    assert_eq!(counts.reference_reads, 2.0);
    assert_eq!(counts.reference_forward_reads, 1.0);
    assert_eq!(counts.reference_placed_left_reads, 1.0);
    // Every allele's reads, both samples: 2 + 4 + 5 in the first, 2 + 4 in the second.
    assert_eq!(counts.total_reads, 17.0);
    // **And the ninth count is read from the primary alternative's copies, not from
    // allele 1's.** Both samples are called homozygous for allele 2, so the calls expect
    // every one of the 17 reads to carry it; counting copies of allele 1 — which neither
    // call carries — would give 0, a maximal apparent deficit at a locus with none.
    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(2), AlleleId(2)]
    );
    assert_eq!(counts.genotype_expected_alternative_reads, 17.0);
}

/// **Two alternatives the same number of reads reached: the one whose bases sort first
/// wins, whatever its id.**
///
/// Ids follow the order the merge met the alleles, which is the order the samples were
/// given, so an id tie-break changed a site's quality with the cohort's order. Here allele
/// 1 is `T` and allele 2 is `C`, so the lower id and the earlier bases disagree, and no
/// sample carries either, so the choice falls back to all alternatives.
#[test]
fn a_tie_between_two_alternatives_goes_to_the_bases_that_sort_first() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    let rows = [
        observation(0, 0, 2, 1, 1),
        observation(1, 0, 5, 3, 2),
        observation(2, 0, 5, 1, 4),
    ];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![0.0, -1.0, -3.0, -1.0, -3.0, -3.0]],
        &evidence,
        &[1.0, 0.5, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(0)]
    );
    let counts = inference
        .artifact_test_counts()
        .expect("both alternatives drew reads");
    assert_eq!(counts.primary_alternative, AlleleId(2), "C sorts before T");
    // And the counts are that allele's, not the other's — 1 forward against 3.
    assert_eq!(counts.alternative_forward_reads, 1.0);
}

/// **The most reads win over a lower id, and a tie is not won by the higher id.**
///
/// The fixture above ties two alternatives whose higher id is also the one whose bases
/// sort first, so a rule that let a later allele win ties, or one that replaced the best
/// only on a tie, would pass it. Here allele 1 is `T` with 2 reads, and alleles 2 (`C`) and
/// 3 (`G`) draw 5 each: the most reads rule out `T` although it comes first, and of the two
/// tied `C` wins although `G` comes after it. No sample carries any of them.
#[test]
fn more_reads_beat_a_lower_id_and_a_tie_goes_to_the_bases_not_the_later_id() {
    let (alleles, table) = generic_locus(3);
    let view = table.view();
    let rows = [
        observation(0, 0, 9, 4, 5),
        observation(1, 0, 2, 1, 1),
        observation(2, 0, 5, 2, 3),
        observation(3, 0, 5, 4, 1),
    ];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let mut favouring_reference = vec![-6.0; view.genotype_count()];
    favouring_reference[0] = 0.0;
    let inference = call_locus(
        &[favouring_reference],
        &evidence,
        &[1.0, 0.5, 0.5, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(0)]
    );
    let counts = inference
        .artifact_test_counts()
        .expect("all three alternatives drew reads");
    assert_eq!(
        counts.primary_alternative,
        AlleleId(2),
        "C: 5 reads, and before G"
    );
    assert_eq!(counts.alternative_forward_reads, 2.0);
}

/// **An alternative a sample is called carrying is weighed before one with more reads that
/// no sample carries.**
///
/// The tests are about the record: the allele-balance test expects the called copies of
/// the alternative it weighs, and an allele absent from every call expects none and is not
/// in the record's `ALT`. On 8 tomato accessions, SL4.0ch01:3493044 lost 126 of 163 Phred
/// to the strand test on such an allele. Here allele 1 draws 5 reads and no call; allele 2
/// draws 2 and the sample is called `0/2`.
#[test]
fn the_primary_alternative_is_one_a_call_carries_before_one_with_more_reads() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    let rows = [
        observation(0, 0, 4, 2, 2),
        observation(1, 0, 5, 5, 0),
        observation(2, 0, 2, 1, 1),
    ];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    // Genotypes in table order 0/0, 0/1, 1/1, 0/2, 1/2, 2/2; the reads favour 0/2.
    let inference = call_locus(
        &[vec![-6.0, -6.0, -6.0, 0.0, -6.0, -6.0]],
        &evidence,
        &[1.0, 0.5, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(2)]
    );
    let counts = inference
        .artifact_test_counts()
        .expect("both alternatives drew reads");
    assert_eq!(counts.primary_alternative, AlleleId(2));
    assert_eq!(counts.alternative_reads, 2.0);
    // Half of the sample's 11 reads, from the one copy of allele 2 its call carries.
    assert_eq!(counts.genotype_expected_alternative_reads, 5.5);
}

/// **A locus whose alternatives no read reached carries no artifact summary**, and the
/// site quality and the calls are unaffected.
///
/// Both artifact tests weigh one alternative against the reference, so a locus with no
/// alternative *reads* leaves them nothing to weigh; production hands back its baseline
/// unchanged in exactly this case. **Not a rare shape**: the merge builds a locus when a
/// sample's non-reference reads pooled reach its rule, so an alternative that ends the
/// pass with no read of its own is an ordinary outcome.
#[test]
fn a_locus_whose_alternatives_no_read_reached_carries_no_artifact_summary() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 7, 4, 3)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![0.0, -4.0, -9.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert!(inference.artifact_test_counts().is_none());
    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(0)]
    );
    // The site quality is still computed and is the right one for a locus this
    // reference-looking: 8 in 100,000 of a Phred, which is a probability of no variant
    // of 0.99998.
    assert!(
        inference.uncorrected_site_quality().get() < 1e-3,
        "one sample whose reads favour the homozygous reference by 4 nats leaves \
         essentially no chance of a variant here, and the quality says so: {}",
        inference.uncorrected_site_quality().get()
    );
}

/// **A locus called over the reference alone is a first-class result** — 27.4% of built
/// loci on the 63-accession tomato panel and 27.3% on HG002 at 30×
/// (`SelectionVerdict::Selected`) — and it carries no artifact summary, because there is
/// no alternative for the two tests to name.
///
/// One allele means one genotype, so every sample's posterior is `1.0` on it: the
/// genotype quality is the clamp's, capped at 99, rather than an infinity.
#[test]
fn a_locus_called_over_the_reference_alone_carries_no_artifact_summary() {
    let alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
    let table = GenotypeTable::build(diploid(), 1);
    let view = table.view();
    assert_eq!(view.genotype_count(), 1);
    let rows = [observation(0, 0, 5, 3, 2)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![0.0]],
        &evidence,
        &[1.0],
        &alleles,
        &view,
        &FixedLogPriors(vec![0.0]),
    );

    assert!(inference.artifact_test_counts().is_none());
    let (genotype, quality) = called(&inference, 0);
    assert_eq!(genotype.alleles(), [AlleleId(0), AlleleId(0)]);
    assert_eq!(quality.get(), MAX_GENOTYPE_QUALITY);
}

/// **A sample the candidate step set aside is emitted missing, with no quality beside
/// it, and is in neither the artifact counts nor the expectation they are weighed
/// against.**
///
/// Its reads are deliberately loud — 20 of them on the alternative, more than the other
/// two samples together — so that counting them would be visible: with the set-aside
/// sample included the alternative would draw **26 reads of 38** rather than the 6 of 12
/// asserted here, and the expectation would move with it.
///
/// **Why the fixture has to write that sample's likelihood row at all.** Nothing sets a
/// sample aside before the loop yet — that is step D1's — so this test hands the loop a
/// complete table and flags the sample only in the evidence. What it pins is the final
/// pass's half of the ruling, which is the half this step owns.
/// **A sample whose reads said nothing is called and counted nowhere**
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.2): adding one leaves the cohort's copies,
/// and every other sample's genotype and quality, exactly as they were. Before the ruling
/// its prior-only copies entered the sum and moved all three — on the GIAB trio at 30x, `AF`
/// at 766 of 767 records and `GQ` by up to 24, with no read behind any of it.
#[test]
fn a_sample_with_no_reads_changes_no_other_samples_call() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let reads = [
        [observation(0, 0, 6, 3, 3), observation(1, 0, 4, 2, 2)],
        [observation(0, 0, 9, 4, 5), observation(1, 0, 1, 1, 0)],
    ];
    let likelihoods = vec![vec![-3.0, 0.0, -4.0], vec![0.0, -1.5, -7.0]];
    let observed = [called_sample(&reads[0]), called_sample(&reads[1])];
    let alone = call_locus(
        &likelihoods,
        &LocusEvidence::generic(locus_region(), &observed),
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );
    let with_an_empty_sample = [
        called_sample(&reads[0]),
        called_sample(&[]),
        called_sample(&reads[1]),
    ];
    let flat = vec![0.0, 0.0, 0.0];
    let together = call_locus(
        &[likelihoods[0].clone(), flat, likelihoods[1].clone()],
        &LocusEvidence::generic(locus_region(), &with_an_empty_sample),
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert_eq!(
        together.cohort_expected_copies().copies(),
        alone.cohort_expected_copies().copies(),
        "the empty sample put copies in the cohort's sum"
    );
    assert_eq!(together.per_sample[0], alone.per_sample[0]);
    assert_eq!(together.per_sample[2], alone.per_sample[1]);
    assert!(matches!(
        together.per_sample[1],
        SampleGenotypeCall::Called {
            reads_were_uninformative: true,
            ..
        }
    ));
}

#[test]
fn a_sample_the_candidate_step_set_aside_is_missing_and_counted_nowhere() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let quiet = [observation(0, 0, 3, 2, 1), observation(1, 0, 3, 1, 2)];
    let loud = [observation(0, 0, 6, 3, 3), observation(1, 0, 20, 18, 2)];
    let per_sample = [
        called_sample(&quiet),
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(&loud, 0.0, &[]),
            genotype_must_be_missing: true,
        },
        called_sample(&quiet),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let heterozygous = vec![-2.0, 0.0, -2.0];
    let inference = call_locus(
        &[heterozygous.clone(), heterozygous.clone(), heterozygous],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert!(inference.per_sample[1].is_missing());
    assert!(
        inference.per_sample[1].genotype().is_none()
            && inference.per_sample[1].score_best_genotype().is_none(),
        "a sample set aside has no genotype and no quality: there is nothing that was \
         scored to have one"
    );
    assert!(!inference.per_sample[0].is_missing() && !inference.per_sample[2].is_missing());

    let counts = inference
        .artifact_test_counts()
        .expect("the two called samples' reads reached the alternative");
    assert_eq!(counts.total_reads, 12.0, "6 reads from each called sample");
    assert_eq!(counts.alternative_reads, 6.0);
    assert_eq!(counts.reference_reads, 6.0);
    // Both called samples are heterozygous, so the calls expect half of each one's six
    // reads to carry the alternative.
    assert_eq!(counts.genotype_expected_alternative_reads, 6.0);
}

/// **The expectation is read from the calls and not from the fitted frequency**, which
/// is what makes the allele-balance test able to catch an artifact at all: the frequency
/// adapts to the artifact and would excuse it (`spec/calling_quality.md` §6.2).
///
/// Two samples with the same depth and different calls: one heterozygous, one homozygous
/// reference. The calls expect `0.5 × 8 + 0 × 8 = 4` alternative reads. The loop's
/// converged copies at this locus are `[3.0441, 0.9559]`, so the fitted alternative
/// frequency is **0.23898** of the cohort's four chromosomes and a version that used it
/// would expect `0.23898 × 16 = 3.8237` — close enough to be mistaken for the right
/// answer and not equal to it, which is the point of asserting the exact 4.
#[test]
fn the_expected_alternative_reads_come_from_the_calls_and_not_from_the_frequency() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let heterozygote = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let homozygous_reference = [observation(0, 0, 7, 4, 3), observation(1, 0, 1, 1, 0)];
    let per_sample = [
        called_sample(&heterozygote),
        called_sample(&homozygous_reference),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![-4.0, 0.0, -4.0], vec![0.0, -4.0, -9.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    assert_eq!(
        called(&inference, 1).0.alleles(),
        [AlleleId(0), AlleleId(0)]
    );
    let counts = inference.artifact_test_counts().expect("a summary");
    assert_eq!(counts.total_reads, 16.0);
    assert_eq!(counts.genotype_expected_alternative_reads, 4.0);
}

/// **At ploidy four a call carrying one copy expects a quarter of its reads**, not a
/// half — the divisor is the locus's ploidy and the fixture is the only one here that
/// can tell the two apart.
///
/// Every other fixture in this module is diploid, where the ploidy, the allele count and
/// the copies of a heterozygote are all 2 and a hard-coded divisor passes. Here the winning genotype is `0/0/0/1`: the minted call names four alleles,
/// and the expectation is `(1 ÷ 4) × 8 = 2`, against the 4 a divisor of two would give.
#[test]
fn at_ploidy_four_a_call_carrying_one_copy_expects_a_quarter_of_its_reads() {
    let mut alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
    alleles.admit(Box::from(b"T".as_slice()));
    let tetraploid = Ploidy::try_new(4).expect("a tetraploid");
    let table = GenotypeTable::build(tetraploid, 2);
    let view = table.view();
    // Five genotypes: 4/0, 3/1, 2/2, 1/3, 0/4 copies of (reference, alternative).
    assert_eq!(view.genotype_count(), 5);
    assert_eq!(view.allele_counts_of(GenotypeIdx(1)), Some(&[3, 1][..]));

    let rows = [observation(0, 0, 6, 3, 3), observation(1, 0, 2, 1, 1)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![-3.0, 0.0, -3.0, -6.0, -9.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &FixedLogPriors(vec![0.0; 5]),
    );

    let (genotype, _) = called(&inference, 0);
    assert_eq!(
        genotype.alleles(),
        [AlleleId(0), AlleleId(0), AlleleId(0), AlleleId(1)],
        "a tetraploid call names one allele per copy of the genome"
    );
    let counts = inference.artifact_test_counts().expect("a summary");
    assert_eq!(counts.total_reads, 8.0);
    assert_eq!(counts.genotype_expected_alternative_reads, 2.0);
}

/// **A read that saw only part of the locus is in neither the depth nor the total**, and
/// the reason is not tidiness: a partial says the sample carries *at least* this, not
/// what it carries, so it stands behind no allele and neither artifact test has a column
/// for it.
///
/// The fixture's partial carries 50 reads against the 8 the summary counts — more than
/// six times as many — so a version that folded partials into the depth would miss by a
/// factor rather than by a rounding.
#[test]
fn reads_that_saw_only_part_of_the_locus_are_in_neither_the_depth_nor_the_total() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
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

    let inference = call_locus(
        &[vec![-4.0, -4.0, 0.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    let counts = inference.artifact_test_counts().expect("a summary");
    assert_eq!(counts.total_reads, 8.0);
    assert_eq!(counts.genotype_expected_alternative_reads, 8.0);
}

/// **The cohort's expected copies leave the locus as the loop left them, not recomputed
/// from the calls** (`spec/calling_em_loop.md` §9).
///
/// Two samples, both called heterozygous. Their calls say the cohort carries two copies
/// of each allele; the loop's converged estimate is **`[2.6759, 1.3241]`**, which is
/// what the record carries. **Two thirds of a copy apart**, and the gap is the
/// uncertainty a call has thrown away: each sample's converged posterior is
/// `[0.3438, 0.6507, 0.0055]`, so **a third** of its probability sits on the homozygous
/// reference — 2 × 0.3438 is the two thirds of a copy — and a call has no way to say it.
/// Site filtering and emission both read that number.
///
/// Asserted **bit for bit** against a second run of the same loop, so the check is that
/// the value travelled rather than that it is approximately right.
#[test]
fn the_cohort_expected_copies_are_the_loops_and_not_recomputed_from_the_calls() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let per_sample = [called_sample(&rows), called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    // Reads that favour the heterozygote by one nat over the homozygous reference — so
    // each sample keeps real probability on a genotype its call does not name.
    let likelihoods = vec![vec![-1.0, 0.0, -4.0], vec![-1.0, 0.0, -4.0]];

    let inference = call_locus(
        &likelihoods,
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    let (_, loop_only) = loop_over(
        &likelihoods,
        &[1.0, 0.5],
        &alleles,
        &view,
        &RunnableCallingLoopConfig::default(),
    );
    assert_eq!(
        inference.cohort_expected_copies().copies(),
        loop_only.cohort_expected_copies(),
        "the record carries the loop's own row"
    );
    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    assert_eq!(
        called(&inference, 1).0.alleles(),
        [AlleleId(0), AlleleId(1)]
    );
    let copies = inference.cohort_expected_copies().copies();
    assert!(
        (copies[0] - 2.675_948_480_830_401).abs() < 1e-12
            && (copies[1] - 1.324_051_519_169_599).abs() < 1e-12,
        "the loop's converged row is [2.6759, 1.3241] and the record carries {copies:?}"
    );
    assert!(
        (copies[0] - 2.0).abs() > 0.5,
        "two heterozygous calls would give exactly two copies of each allele, and the \
         expected copies are two thirds of a copy away from that — which is the \
         uncertainty the calls threw away"
    );
}

/// **The site quality on the record is the fold's own number for this likelihood
/// table** — the same value a second, independently prepared scratch produces from the
/// same rows.
///
/// What this catches is the field arriving unwritten or written from something else: a
/// quality of zero is a perfectly ordinary answer at a locus nobody carries, so a test
/// that only asked for *a* number would pass against a pass that lost it. This fixture's
/// two samples both favour the homozygous alternative by 9 nats, and the quality is
/// **42.373**.
#[test]
fn the_site_quality_on_the_record_is_the_fold_of_the_table_the_loop_built() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 1), observation(1, 0, 6, 5, 2)];
    let per_sample = [called_sample(&rows), called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let likelihoods = vec![vec![-9.0, -4.0, 0.0], vec![-9.0, -4.0, 0.0]];

    let inference = call_locus(
        &likelihoods,
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    let mut fold_only = CallingScratch::<()>::default();
    fold_only.prepare_for_locus(likelihoods.len(), &alleles, &view);
    for (sample, row) in likelihoods.iter().enumerate() {
        for (slot, &value) in fold_only
            .sample_genotype_likelihoods_mut(sample)
            .iter_mut()
            .zip(row)
        {
            *slot = LogProb(value);
        }
    }
    let expected = score_uncorrected_site_quality(
        fold_only.site_quality_buffers_mut(),
        &view,
        human_like_seed(),
    );
    assert_eq!(inference.uncorrected_site_quality(), expected);
    assert!(
        expected.get() > 0.0,
        "two samples whose reads exclude the reference give a quality well above zero, \
         so this fixture can tell a lost value from a computed one"
    );
}

/// **A repeat tract carries a site quality and no artifact summary**, and no sample there
/// is set aside.
///
/// The strand and read-position tests are the SNP/indel path's: at a tract what goes
/// wrong is slippage, which is already inside the read likelihood, and
/// `spec/calling_quality.md` §8 leaves a tract's quality to a sibling document that is
/// not written. **The likelihood table here is the fixture's own**: the repeat-tract row
/// exists (`likelihood/ssr.rs`), but the repeat-tract *candidate* path does not, so
/// nothing upstream of this pass yet hands a tract the alleles to score.
#[test]
fn a_repeat_tract_carries_a_site_quality_and_no_artifact_summary() {
    let detail = SsrDetail {
        motif: Motif::new(b"AT").expect("a dinucleotide motif"),
        left_flank: Box::from(b"CCCGGG".as_slice()),
        right_flank: Box::from(b"TTTAAA".as_slice()),
    };
    let mut alleles = CandidateAlleles::new(
        Box::from(b"ATAT".as_slice()),
        LocusKind::Ssr(SsrDetail {
            motif: Motif::new(b"AT").expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        }),
    );
    alleles.admit(Box::from(b"ATATAT".as_slice()));
    let table = GenotypeTable::build(diploid(), 2);
    let view = table.view();
    let per_sample = [SsrSampleEvidence::new(&[], &detail, TEST_TRACT_LEN)];
    let repeat_counts = tract_repeat_counts();
    let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);

    let inference = call_locus(
        &[vec![-4.0, 0.0, -4.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    assert!(inference.artifact_test_counts().is_none());
    assert!(!inference.per_sample[0].is_missing());
    assert_eq!(
        called(&inference, 0).0.alleles(),
        [AlleleId(0), AlleleId(1)]
    );

    // **The site quality is the only number this arm has to get right**, so it is
    // compared against an independently prepared fold rather than merely being present —
    // and the fold runs under the tract's own seed (1.0 for the reference tract against
    // 0.5 for the other), not the run's (`calling_quality.md` §8.1).
    let mut fold_only = CallingScratch::<()>::default();
    fold_only.prepare_for_locus(1, &alleles, &view);
    for (slot, &value) in fold_only
        .sample_genotype_likelihoods_mut(0)
        .iter_mut()
        .zip(&[-4.0_f64, 0.0, -4.0])
    {
        *slot = LogProb(value);
    }
    let expected = score_uncorrected_site_quality(
        fold_only.site_quality_buffers_mut(),
        &view,
        SiteQualityPrior::of_repeat_tract(&[1.0, 0.5]),
    );
    let under_the_run_seed = score_uncorrected_site_quality(
        fold_only.site_quality_buffers_mut(),
        &view,
        human_like_seed(),
    );
    assert!(expected.get() > 0.0, "the fixture's fold is above zero");
    assert!(
        expected.get() > under_the_run_seed.get() + 10.0,
        "a tract seed expecting variation must ask for less evidence than a seed of 1 \
         variant in 1,000: {} against {}",
        expected.get(),
        under_the_run_seed.get()
    );
    assert_eq!(inference.uncorrected_site_quality(), expected);
}

/// Evidence and the run's inbreeding coefficients covering different cohorts are refused:
/// the pass walks the coefficients and indexes the evidence by the same number, so two
/// lists of different lengths are two different sample orders.
#[test]
#[should_panic(expected = "one per sample of the run")]
fn evidence_covering_another_cohort_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 1)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = FrozenParameters::uncontaminated(
        &calibration,
        &inbreeding,
        human_like_seed(),
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &view);
    let _ = summarise_final_pass(
        &mut scratch,
        &view,
        &evidence,
        &parameters,
        &MarginalizedDirichletPrior,
        alleles.clone(),
        FrequencyLoopOutcome {
            passes: 1,
            converged: true,
        },
        Provenance::FittedHere,
        None,
    );
}

/// One inbreeding coefficient per sample of the run, **in both directions** — the walk is
/// over the coefficients, so a short slice would call some samples and leave the rest with
/// no call at all.
#[test]
#[should_panic(expected = "one per sample of the run")]
fn fewer_inbreeding_coefficients_than_samples_are_refused_by_the_final_pass() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 1)];
    let per_sample = [called_sample(&rows), called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = FrozenParameters::uncontaminated(
        &calibration,
        &inbreeding,
        human_like_seed(),
        &strata,
        &NO_SUBSTITUTION_RATES,
        diploid(),
    );
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &view);
    let _ = summarise_final_pass(
        &mut scratch,
        &view,
        &evidence,
        &parameters,
        &MarginalizedDirichletPrior,
        alleles.clone(),
        FrequencyLoopOutcome {
            passes: 1,
            converged: true,
        },
        Provenance::FittedHere,
        None,
    );
}

/// **The callable samples and the scratch's rows must be the same set**, and the row map
/// is what says so. A scratch prepared for two rows at a locus with one callable sample
/// was filled for a different set of them — and every row it did read was a legal row of a
/// real table, so nothing about the arithmetic complains. The join is checked at each
/// sample rather than by a count at the end, because a count is satisfied by a
/// *permutation*: the table filled for one sample and read for another.
#[test]
#[should_panic(expected = "was filled for a different set of samples")]
fn a_scratch_prepared_for_more_rows_than_the_locus_can_call_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 1), observation(1, 0, 6, 5, 2)];
    let per_sample = [
        called_sample(&rows),
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(&rows, 0.0, &[]),
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
        view.ploidy(),
    );
    let mut scratch = CallingScratch::<()>::default();
    // Two rows where the locus can call one, and both claimed — the shape the driver
    // cannot produce and a hand-assembled call site can.
    scratch.prepare_for_locus(2, &alleles, &view);
    scratch.claim_row_for(0, outbred());
    scratch.claim_row_for(1, outbred());
    for row in 0..2 {
        for (slot, &value) in scratch
            .sample_genotype_likelihoods_mut(row)
            .iter_mut()
            .zip(&[-4.0_f64, 0.0, -4.0])
        {
            *slot = LogProb(value);
        }
    }
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[2.0, 2.0]);
    for row in 0..2 {
        scratch
            .sample_expected_copies_mut(row)
            .copy_from_slice(&[1.0, 1.0]);
    }
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    let _ = summarise_final_pass(
        &mut scratch,
        &view,
        &evidence,
        &parameters,
        &MarginalizedDirichletPrior,
        alleles.clone(),
        FrequencyLoopOutcome {
            passes: 1,
            converged: true,
        },
        Provenance::FittedHere,
        None,
    );
}

/// **Rows claimed past the end of the run are the one disagreement the per-sample check
/// cannot see**, because a row nobody walks is never asked whose it is. Three rows claimed
/// at a two-sample run: every sample matches its own row and the walk still ends one row
/// short of the table, which is the only trace there is.
#[test]
#[should_panic(expected = "callable samples and the scratch was prepared for")]
fn rows_claimed_past_the_end_of_the_run_are_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 1), observation(1, 0, 6, 5, 2)];
    let per_sample = [called_sample(&rows), called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let calibration = [ReadGroupCalibration::unscaled()];
    let inbreeding = [outbred(), outbred()];
    let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
    let parameters = uncontaminated_run(
        &calibration,
        &inbreeding,
        &strata,
        &NO_SUBSTITUTION_RATES,
        view.ploidy(),
    );
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(3, &alleles, &view);
    for run_sample in 0..3 {
        scratch.claim_row_for(run_sample, outbred());
    }
    for row in 0..3 {
        for (slot, &value) in scratch
            .sample_genotype_likelihoods_mut(row)
            .iter_mut()
            .zip(&[-4.0_f64, 0.0, -4.0])
        {
            *slot = LogProb(value);
        }
        scratch
            .sample_expected_copies_mut(row)
            .copy_from_slice(&[1.0, 1.0]);
    }
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[3.0, 3.0]);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    let _ = summarise_final_pass(
        &mut scratch,
        &view,
        &evidence,
        &parameters,
        &MarginalizedDirichletPrior,
        alleles.clone(),
        FrequencyLoopOutcome {
            passes: 1,
            converged: true,
        },
        Provenance::FittedHere,
        None,
    );
}

/// An observation naming an allele the locus was not called over is refused **by name**:
/// without the check the buffer index panics with a message about a slice, which sends
/// the reader to the scratch rather than to the mapping that produced the id.
#[test]
#[should_panic(expected = "mapped against a different allele table")]
fn an_observation_naming_an_allele_the_locus_lacks_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 1), observation(4, 0, 3, 1, 1)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let _ = call_locus(
        &[vec![0.0, -1.0, -3.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );
}

/// More forward-strand reads than reads is refused. Both artifact tests read the two
/// counts as fractions of the read count, and a fraction above one reaches the binomial
/// tail as a probability above one — which comes back as a penalty rather than as a
/// failure.
#[test]
#[should_panic(expected = "shares of those same reads")]
fn more_forward_reads_than_reads_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 3, 1)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let _ = call_locus(
        &[vec![0.0, -1.0, -3.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );
}

/// The placed-left half of the same check, which a fixture varying only the strand count
/// would leave unreached.
#[test]
#[should_panic(expected = "shares of those same reads")]
fn more_placed_left_reads_than_reads_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 5)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let _ = call_locus(
        &[vec![0.0, -1.0, -3.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );
}

/// **Each sample's confidence is its own**, taken from the posterior row while that sample
/// is the one in it — which is why this is a pass and not a read-back after the loop
/// (`spec/calling_quality.md` §3.1). `CallingScratch`'s posterior row is one buffer that
/// every sample in turn is scored into, so a pass that took the quality *after* the walk
/// would give every sample the last one's number.
///
/// Two samples, one favouring the heterozygote by 4 nats and one the homozygous reference
/// by 4, so their winning posteriors are far apart: **11.553 and 17.697 Phred**, six
/// Phred and a call in fifty apart. A pass that gave both samples one sample's quality
/// leaves every other test in this module green.
#[test]
fn two_samples_take_their_confidence_from_their_own_posterior_rows() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 4, 2, 2), observation(1, 0, 4, 2, 2)];
    let per_sample = [called_sample(&rows), called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![-4.0, 0.0, -4.0], vec![0.0, -4.0, -9.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    let first = called(&inference, 0).1.get();
    let second = called(&inference, 1).1.get();
    assert!(
        (first - 11.552_688).abs() < 1e-3 && (second - 17.697_23).abs() < 1e-3,
        "each sample's quality comes from its own posterior row: got {first} and {second}"
    );
}

/// **A set-aside sample's reads do not choose the alternative either**, not merely the
/// counts (`spec/calling_quality.md` §6.3).
///
/// Over the two called samples allele 1 draws 10 reads and allele 2 draws 6; the set-aside
/// sample's 20 reads of allele 2 would reverse that if they were pooled. **The biallelic
/// fixture beside this one cannot catch it** — there a set-aside sample's reads can move
/// the counts but not which allele is counted, because there is only one alternative to
/// choose.
#[test]
fn the_primary_alternative_ignores_a_set_aside_samples_reads() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    let called_rows = [
        observation(0, 0, 2, 1, 1),
        observation(1, 0, 5, 3, 2),
        observation(2, 0, 3, 2, 1),
    ];
    let set_aside_rows = [observation(2, 0, 20, 10, 10)];
    let per_sample = [
        called_sample(&called_rows),
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(&set_aside_rows, 0.0, &[]),
            genotype_must_be_missing: true,
        },
        called_sample(&called_rows),
    ];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let row = vec![0.0, -1.0, -3.0, -1.0, -3.0, -3.0];

    let inference = call_locus(
        &[row.clone(), row.clone(), row],
        &evidence,
        &[1.0, 0.5, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );

    let counts = inference.artifact_test_counts().expect("a summary");
    assert_eq!(
        counts.primary_alternative,
        AlleleId(1),
        "over the called samples allele 1 drew 10 reads and allele 2 drew 6; the \
         set-aside sample's 20 reads of allele 2 are not pooled"
    );
}

/// **The four values the pass only carries** reach the record unchanged: where the locus
/// is, whether the loop converged, how many passes it took, and where its parameters came
/// from. None of the four panics when wrong — a transposed region is a silently misplaced
/// variant, and a hard locus reported as converged is a stronger claim than the loop
/// earned (`spec/calling_em_loop.md` §6).
#[test]
fn the_pass_carries_the_region_the_outcome_and_the_provenance_onto_the_record() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let rows = [observation(0, 0, 2, 1, 1), observation(1, 0, 6, 5, 2)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);
    let likelihoods = vec![vec![-4.0, 0.0, -4.0]];

    let inference = call_locus(
        &likelihoods,
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &MarginalizedDirichletPrior,
    );
    let (outcome, _) = loop_over(
        &likelihoods,
        &[1.0, 0.5],
        &alleles,
        &view,
        &RunnableCallingLoopConfig::default(),
    );

    assert_eq!(inference.region, locus_region());
    assert_eq!(inference.weakest_provenance, Provenance::FittedHere);
    assert_eq!(inference.passes, outcome.passes);
    assert_eq!(inference.converged, outcome.converged);
    assert!(
        outcome.passes > 1,
        "the fixture takes more than one pass, so `passes` is a real number here and not \
         the constant 1"
    );
}

/// **At ploidy one a call is all or nothing** — one allele named, and either every one of
/// the sample's reads expected to carry the alternative or none of them.
///
/// The low end of the ploidy range this caller commits to, and the divisor's other
/// boundary: the module's other fixtures are diploid and tetraploid, so a mint or a
/// divisor that assumes at least two copies passes them all.
#[test]
fn at_ploidy_one_a_call_expects_all_or_none_of_its_reads() {
    let mut alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
    alleles.admit(Box::from(b"T".as_slice()));
    let haploid = Ploidy::try_new(1).expect("a haploid");
    let table = GenotypeTable::build(haploid, 2);
    let view = table.view();
    // Two genotypes at ploidy one: the reference allele, or the alternative.
    assert_eq!(view.genotype_count(), 2);

    let rows = [observation(0, 0, 2, 1, 1), observation(1, 0, 6, 5, 2)];
    let per_sample = [called_sample(&rows)];
    let evidence = LocusEvidence::generic(locus_region(), &per_sample);

    let inference = call_locus(
        &[vec![-4.0, 0.0]],
        &evidence,
        &[1.0, 0.5],
        &alleles,
        &view,
        &FixedLogPriors(vec![0.0, 0.0]),
    );

    let (genotype, _) = called(&inference, 0);
    assert_eq!(
        genotype.alleles(),
        [AlleleId(1)],
        "a haploid call names one allele, not two"
    );
    let counts = inference.artifact_test_counts().expect("a summary");
    assert_eq!(counts.total_reads, 8.0);
    assert_eq!(
        counts.genotype_expected_alternative_reads, 8.0,
        "the one copy is the alternative, so every read is expected to carry it"
    );
}

proptest! {
    /// **The minted genotype is the copy-count row read as a multiset**: allele `a`
    /// appears exactly as many times as the row says, in ascending id order, and the
    /// whole names one allele per copy of the genome.
    ///
    /// A property rather than a case, because the two end-to-end shapes the module's
    /// other fixtures reach — `[0, 2]` at diploid and `[3, 1]` at tetraploid — cannot
    /// separate a mint that drops a zero-copy allele's *position* from one that does not.
    #[test]
    fn mint_genotype_repeats_each_allele_as_many_times_as_the_row_says(
        copies in proptest::collection::vec(0_u32..4, 1..6)
            .prop_filter("a genotype has at least one copy", |row: &Vec<u32>| {
                row.iter().sum::<u32>() > 0
            })
    ) {
        let minted = mint_genotype(&copies);
        let total: u32 = copies.iter().sum();
        proptest::prop_assert_eq!(minted.alleles().len(), total as usize);
        for (allele, &want) in copies.iter().enumerate() {
            let got = minted
                .alleles()
                .iter()
                .filter(|id| usize::from(id.get()) == allele)
                .count();
            proptest::prop_assert_eq!(got, want as usize);
        }
        let mut sorted = minted.alleles().to_vec();
        sorted.sort_by_key(|id| id.get());
        proptest::prop_assert_eq!(&sorted[..], minted.alleles());
    }
}
