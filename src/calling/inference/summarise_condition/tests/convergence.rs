use super::*;

// ---------------------------------------------------------------------------------
// Convergence, the cap, and the emitted flag
// ---------------------------------------------------------------------------------

/// A configuration a locus can meet **only by reaching a bitwise fixed point**, so that
/// short of one the cap is what stops it.
///
/// **Not "never settles", which is what this was first called and is not true of any
/// fixture here.** A threshold of `1e-300` is met the moment two passes produce identical
/// copies, and expectation-maximization reaches that: measured,
/// `three_samples_pulling_apart` hits its fixed point at pass 29 and reports
/// `converged: true` however small the threshold. The helper is honest below that and a
/// trap above it.
///
/// **`1e-300`, not zero** — [`CallingLoopConfig::validate`] refuses a threshold that is
/// not above zero, and it is right to.
fn settles_only_at_a_bitwise_fixed_point(max_passes: u32) -> RunnableCallingLoopConfig {
    CallingLoopConfig {
        max_passes: NonZeroU32::new(max_passes).expect("a cap of at least one pass"),
        convergence_threshold: 1e-300,
        ..CallingLoopConfig::DEFAULT
    }
    .validate()
    .expect("a positive finite threshold below the ceiling")
}

/// Three samples whose reads pull them apart, at a biallelic locus — the fixture the
/// cap tests and the hand-driven oracle share.
///
/// The rows are the three diploid genotypes in the table's order (`0/0`, `0/1`, `1/1`).
/// One sample's reads favour each, by 2 nats — about 8.7 Phred, a handful of reads —
/// which is enough disagreement that the cohort's frequencies are still moving after
/// several passes.
fn three_samples_pulling_apart() -> Vec<Vec<f64>> {
    vec![
        vec![0.0, -2.0, -6.0],
        vec![-2.0, 0.0, -2.0],
        vec![-6.0, -2.0, 0.0],
    ]
}

/// **One sample reaches its fixed point on the loop's first pass, and the second only
/// confirms it** — `spec/calling_em_loop.md` §13 test 1.
///
/// Two runs of the same locus, differing only in the pass cap.
///
/// - **Capped at one pass**, the loop stops with `converged = false` and the two cohort
///   rows it leaves behind are the prior-free initialisation's and the first seeded
///   pass's. §3 says they must differ, because the seed is worth 20 to 30 Phred, and
///   they do: measured on this fixture, `[1.7279, 0.2721]` against `[1.9998, 0.0002]`.
///   That is 0.272 copies of the alternative allele, against the 0.002 raw copies the
///   `1e-3` threshold allows over one diploid sample's two chromosomes — 136 times the
///   movement that would have stopped the loop.
/// - **At the shipped cap** the loop stops after **two** passes, and pass 2's copies are
///   pass 1's **bit for bit** — asserted with `assert_eq!` on the slices, not a
///   tolerance.
///
/// **No branch on cohort size is what produces this** (§7): the prior's leave-one-out
/// term subtracts the one sample's copies from a cohort total that *is* those copies, so
/// the concentration comes back as the seed on every pass and the second pass cannot
/// move. Whether to skip that second pass is spec §12's question 6; this test asserts on
/// the genotype and on pass-1-equals-pass-2 as well, so it survives the branch being
/// added.
#[test]
fn at_one_sample_the_loop_settles_on_the_second_pass() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let seed = [1.0, 0.000_5];
    // Reads that mildly favour the reference homozygote: 1 nat over the heterozygote.
    let likelihoods = vec![vec![0.0, -1.0, -6.0]];

    let (capped, scratch) = loop_over(&likelihoods, &seed, &alleles, &view, &capped_at(1));
    assert_eq!(capped.passes, 1);
    assert!(
        !capped.converged,
        "the prior-free initialisation and the first seeded pass are 20 to 30 Phred \
         apart, so one pass cannot have settled"
    );
    let initialisation = scratch.previous_cohort_expected_copies().to_vec();
    let after_one_pass = scratch.cohort_expected_copies().to_vec();
    assert!(
        (initialisation[1] - after_one_pass[1]).abs() > 0.1,
        "the seed moved the alternative allele's copies from {} to {}, which is less \
         than a tenth of a copy — the fixture no longer poses the question",
        initialisation[1],
        after_one_pass[1]
    );

    let (settled, scratch) = loop_over(
        &likelihoods,
        &seed,
        &alleles,
        &view,
        &RunnableCallingLoopConfig::default(),
    );
    assert!(settled.converged);
    assert_eq!(
        settled.passes, 2,
        "pass 1 differs from the initialisation and pass 2 cannot differ from pass 1"
    );
    assert_eq!(
        scratch.previous_cohort_expected_copies(),
        scratch.cohort_expected_copies(),
        "at one sample pass 2's expected copies are pass 1's bit for bit"
    );
    let (winner, _) = scratch
        .posterior_row()
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.partial_cmp(b).expect("a finite posterior"))
        .expect("a locus has at least one genotype");
    assert_eq!(winner, 0, "the reads and the seed both favour 0/0");
}

/// **A locus that runs out of passes is emitted with the flag set** —
/// `spec/calling_em_loop.md` §13 test 4 and §6.
///
/// The same three-sample locus twice, differing only in the pass cap. Capped at two it
/// reports `passes = 2, converged = false`; at the shipped cap of 50 it settles in 4
/// passes and reports `converged = true`. **So the flag is about the cap and not about
/// the locus**, which is the property that makes it worth emitting.
///
/// **The capped locus is called, and the flag survives the call**: the outcome's two
/// fields go into a `LocusInference` and come back out. `LocusInference::new`
/// deliberately does not check `converged` against `passes` — the cap is run
/// configuration it cannot see — so this is what says a capped locus is emitted rather
/// than refused.
///
/// **Deliberately not asserted: that the delta falls on every pass.**
/// Expectation-maximization guarantees a monotone likelihood, not a monotone parameter
/// delta, and §6 claims no such thing (§13 test 4).
#[test]
fn a_locus_that_runs_out_of_passes_is_called_and_says_so() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let seed = [1.0, 0.000_5];
    let likelihoods = three_samples_pulling_apart();

    let (capped, _) = loop_over(&likelihoods, &seed, &alleles, &view, &capped_at(2));
    assert_eq!(
        capped,
        FrequencyLoopOutcome {
            passes: 2,
            converged: false,
        },
        "two passes is not enough for this locus, so the cap is what stopped it"
    );

    let (settled, _) = loop_over(
        &likelihoods,
        &seed,
        &alleles,
        &view,
        &RunnableCallingLoopConfig::default(),
    );
    assert!(
        settled.converged,
        "the same locus settles when it is allowed the passes"
    );
    assert_eq!(settled.passes, 4);

    // The flag reaches the output, which is the half §6 says nothing downstream can
    // reconstruct.
    let called = LocusInference::new(
        GenomeRegion {
            contig: ContigId(3),
            start: Position(940),
            end: Position(940),
        },
        alleles.clone(),
        vec![
            SampleGenotypeCall::Called {
                genotype: Genotype::new(vec![AlleleId(0), AlleleId(0)]),
                genotype_quality: Phred::try_new(30.0).expect("a quality"),
                reads_were_uninformative: false,
                genotype_rows: crate::calling::placeholder_genotype_rows(),
            };
            3
        ],
        ExpectedAlleleCopies::new(vec![5.0, 1.0], &alleles),
        capped.converged,
        capped.passes,
        Provenance::FittedHere,
        None,
        a_worker_written_site_quality(),
        None,
    );
    assert!(!called.converged);
    assert_eq!(called.passes, 2);
}

/// **The same movement on the frequency scale stops the loop at one sample and at a
/// thousand** — the division's own test (`spec/calling_em_loop.md` §6).
///
/// Expected copies are a count and the threshold is a fraction, so the criterion is
/// written on the count *divided by the cohort's chromosomes*. This walks a cohort of 1
/// and a cohort of 1,000, gives each the raw movement that is the same fraction of its
/// own chromosome total, and requires the same verdict: settled at 9 parts in 10,000,
/// not settled at 11.
///
/// **Measured against the mutation it exists for.** With the `/ cohort_chromosomes`
/// deleted, this test stops at its very first cell: at one sample the movement is 0.0018
/// raw copies against a `1e-3` threshold, so a locus that *has* settled is reported as
/// still moving. The thousand-sample cell, at 1.8 raw copies, would fail the same way —
/// which is the point. **A criterion on raw counts is not merely tighter at a thousand
/// samples; it is already wrong at one.** The same mutation also pushes the cap test's
/// three-sample locus from 4 passes to 6.
#[test]
fn the_same_frequency_scale_movement_settles_at_one_sample_and_at_a_thousand() {
    let threshold = DEFAULT_CONVERGENCE_THRESHOLD;
    for samples in [1_usize, 1000] {
        let chromosomes = 2.0 * samples as f64;
        for (frequency_movement, should_settle) in [(9e-4, true), (1.1e-3, false)] {
            let previous = [10.0, 0.0];
            let current = [10.0, frequency_movement * chromosomes];
            assert_eq!(
                cohort_expected_copies_have_settled(
                    &previous,
                    &current,
                    diploid(),
                    samples,
                    threshold
                ),
                should_settle,
                "a movement of {frequency_movement} on the frequency scale is \
                 {} raw copies over {chromosomes} chromosomes",
                current[1]
            );
        }
    }
}

/// **Three of the four natural spellings of the delta report a locus settled before any
/// pass has advanced**, and the shipped one does not.
///
/// `prepare_for_locus` fills **both** cohort rows with `UNWRITTEN_SCRATCH_VALUE`, which
/// is `NaN`, and the previous row keeps it until a pass has actually swapped. Every
/// comparison against a `NaN` is false, and `f64::max` is documented to return *the other
/// argument* when one is `NaN` — so a fold discards the sentinel instead of propagating
/// it, and the three fold spellings hand back a delta below any threshold.
///
/// The failure that would produce is silent and it is the worst kind: a genotype,
/// **flagged as converged**, from a loop that ran one pass and compared it against
/// nothing.
///
/// This test computes all four on that exact state and asserts the disagreement, so the
/// reason for `all(… < threshold)` is a fact in the suite rather than a remark in a doc
/// comment.
#[test]
fn the_fold_spellings_of_the_delta_settle_where_this_one_does_not() {
    let previous = [UNWRITTEN_SCRATCH_VALUE, UNWRITTEN_SCRATCH_VALUE];
    let current = [1.9, 0.1];
    let chromosomes = 6.0;
    let threshold = DEFAULT_CONVERGENCE_THRESHOLD;
    let scaled: Vec<f64> = previous
        .iter()
        .zip(&current)
        .map(|(before, after)| (after - before).abs() / chromosomes)
        .collect();

    let from_zero = scaled.iter().copied().fold(0.0_f64, f64::max);
    let from_negative_infinity = scaled.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut by_hand = 0.0_f64;
    for &movement in &scaled {
        if movement > by_hand {
            by_hand = movement;
        }
    }

    assert!(
        from_zero < threshold,
        "fold(0.0, f64::max) discards the sentinel and returns {from_zero}"
    );
    assert!(
        from_negative_infinity < threshold,
        "fold(-inf, f64::max) discards it too and returns {from_negative_infinity}"
    );
    assert!(
        by_hand < threshold,
        "a hand-written `>` comparison loses every comparison against a NaN and \
         returns {by_hand}"
    );
    assert!(
        !cohort_expected_copies_have_settled(&previous, &current, diploid(), 3, threshold),
        "the shipped rule is the only one of the four that refuses to settle against a \
         row no pass has written"
    );
}

/// **The loop is exactly the sequence of passes driven by hand**, cohort row for cohort
/// row, bit for bit.
///
/// The oracle for everything the loop's own body decides and no single-property test
/// pins: that the prior-free pass runs once and is not counted, that the E-step is
/// scored against the row the *previous* M-step wrote, that the swap sits between the
/// two halves of a pass, and that `passes` counts the seeded passes.
///
/// **The two sides are stopped by different means and that is deliberate.** The
/// hand-driven side has no stopping rule at all — it is `for _ in 0..PASSES` — so it
/// cannot inherit a bug from the rule under test; the loop is held to the same pass count
/// by a threshold small enough that only a bitwise fixed point meets it, which this
/// fixture does not reach until pass 29. Two rows are compared, not one: the loop's
/// current row against the hand-driven row after four passes, and the loop's *previous*
/// row against the hand-driven row after three — which is what says the swap happens
/// where it does rather than a pass late.
///
/// **Measured against the reordering it exists for.** Moving the swap ahead of the
/// E-step scores every sample against a cohort row that has just been refilled with the
/// `NaN` sentinel. In debug that panics inside the prior — *"the cohort's expected
/// allele copies … got `[NaN, NaN]`"* — but that check is a `debug_assert!`, so
/// **under `--release` the run does not panic**: the leave-one-out `max(0, ·)` returns
/// the other operand on a `NaN`, every sample is scored against the bare seed, and the
/// loop reports `passes: 2, converged: true` where the shipped order gives
/// `passes: 4, converged: false`. A converged flag, two passes early, with the cohort's
/// evidence silently absent. This test is what sees it, and it sees it in both profiles.
#[test]
fn the_loop_reproduces_the_passes_driven_by_hand() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let seed = [1.0, 0.000_5];
    let likelihoods = three_samples_pulling_apart();
    const PASSES: u32 = 4;

    let (outcome, scratch) = loop_over(
        &likelihoods,
        &seed,
        &alleles,
        &view,
        &settles_only_at_a_bitwise_fixed_point(PASSES),
    );
    assert_eq!(
        outcome,
        FrequencyLoopOutcome {
            passes: PASSES,
            converged: false,
        }
    );

    // The same passes, written out.
    let mut by_hand = CallingScratch::<()>::default();
    by_hand.prepare_for_locus(likelihoods.len(), &alleles, &view);
    by_hand.seed_concentration_mut().copy_from_slice(&seed);
    for (sample, row) in likelihoods.iter().enumerate() {
        for (slot, &value) in by_hand
            .sample_genotype_likelihoods_mut(sample)
            .iter_mut()
            .zip(row)
        {
            *slot = LogProb(value);
        }
    }
    for sample in 0..likelihoods.len() {
        score_one_sample(
            by_hand.sample_scoring_buffers_mut(sample),
            &view,
            PassPrior::Flat,
        );
    }
    sum_cohort_expected_copies(by_hand.cohort_summing_buffers_mut());
    let mut rows = Vec::with_capacity(PASSES as usize);
    for _ in 0..PASSES {
        for sample in 0..likelihoods.len() {
            score_one_sample(
                by_hand.sample_scoring_buffers_mut(sample),
                &view,
                PassPrior::LeaveOneOut {
                    model: &MarginalizedDirichletPrior,
                    inbreeding: outbred(),
                },
            );
        }
        by_hand.advance_cohort_expected_copies();
        sum_cohort_expected_copies(by_hand.cohort_summing_buffers_mut());
        rows.push(by_hand.cohort_expected_copies().to_vec());
    }

    assert_eq!(
        scratch.cohort_expected_copies(),
        &rows[PASSES as usize - 1][..],
        "the loop's final row is not the hand-driven pass {PASSES} row"
    );
    assert_eq!(
        scratch.previous_cohort_expected_copies(),
        &rows[PASSES as usize - 2][..],
        "the loop's previous row is not the hand-driven pass {} row, so the swap is \
         not where it should be",
        PASSES - 1
    );
}

/// **Each sample is scored against its own inbreeding coefficient, not the cohort's
/// first** — the property the length check exists to protect, and the one no other
/// fixture reaches, because every other cohort here is `vec![outbred(); n]`.
///
/// Samples 1 and 2 pull toward opposite homozygotes, so trading their coefficients
/// **without moving the samples** must move the answer; an implementation that hands
/// every sample one shared coefficient returns the same row for both runs. Measured, with
/// the coefficient hoisted to `inbreeding_by_sample[0]`: `passes` goes 3 → 4 and the
/// cohort's copies of the alternative allele move 0.38 out of six chromosomes, an
/// allele-frequency shift of 0.064 — finite, plausible, wrong, and until this test,
/// unseen.
///
/// **Not asserted: that moving a sample and its coefficient together reproduces the
/// cohort row bitwise.** It does not, and it should not — the M-step sums in sample order
/// (spec §8), so a permutation moves the sum by one unit in the last place, measured
/// `2.883168867243733` against `2.8831688672437323`. Spec §13 test 2 pins the bitwise
/// comparison on the fixed order for exactly this reason.
#[test]
fn each_sample_is_scored_against_its_own_inbreeding_coefficient() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let seed_concentration = [1.0, 0.000_5];
    let likelihoods = three_samples_pulling_apart();
    let half = InbreedingF::try_new(0.5).expect("a coefficient");
    let most = InbreedingF::try_new(0.9).expect("a coefficient");
    let shipped = RunnableCallingLoopConfig::default();

    let (_, as_given) = loop_over_inbred(
        &likelihoods,
        &[outbred(), half, most],
        &seed_concentration,
        &alleles,
        &view,
        &shipped,
    );
    let (_, traded) = loop_over_inbred(
        &likelihoods,
        &[outbred(), most, half],
        &seed_concentration,
        &alleles,
        &view,
        &shipped,
    );

    let moved = (as_given.cohort_expected_copies()[1] - traded.cohort_expected_copies()[1]).abs();
    assert!(
        moved > 0.01,
        "samples 1 and 2 pull opposite ways, so trading their coefficients must move the \
         cohort's copies of the alternative allele — it moved {moved}, and an \
         implementation that hands every sample one shared coefficient would move it by \
         nothing at all"
    );
}

/// **A locus that settles on its last allowed pass is converged, not capped** — the one
/// input that separates the loop's two exit tests.
///
/// Everywhere else in this file the cap and the settling are several passes apart, so
/// nothing distinguishes the two orderings. `three_samples_pulling_apart` settles on pass
/// 4, so a cap of exactly 4 is where they disagree: measured, testing the cap first
/// returns `converged: false` for a locus that did settle, which understates every
/// genotype at the site (§6).
#[test]
fn run_frequency_loop_reports_converged_when_the_last_allowed_pass_settles() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let (outcome, _) = loop_over(
        &three_samples_pulling_apart(),
        &[1.0, 0.000_5],
        &alleles,
        &view,
        &capped_at(4),
    );
    assert_eq!(
        outcome,
        FrequencyLoopOutcome {
            passes: 4,
            converged: true,
        },
        "this locus settles on pass 4, so a cap of exactly 4 must report it settled \
         rather than capped"
    );
}

/// **The `.abs()` earns its place only from three alleles on**, and this is the row that
/// shows it.
///
/// At a biallelic locus the expected copies sum to the cohort's chromosome total, so the
/// two movements are equal and opposite and a signed comparison gives the same verdict as
/// an absolute one — which is why every other fixture in this file, all of them two
/// alleles, leaves a dropped `.abs()` invisible. At three alleles the reference can fall
/// by more than the threshold while both alternatives rise by less: scaled,
/// `[−0.0015, +0.0008, +0.0007]` against `1e-3`. That is the ordinary shape of a
/// multi-allelic locus still moving, and a caller without the `.abs()` calls it settled.
#[test]
fn a_fall_larger_than_every_rise_has_not_settled() {
    let previous = [10.0, 2.0, 3.0];
    let current = [10.0 - 0.009, 2.0 + 0.004_8, 3.0 + 0.004_2];
    assert!(
        !cohort_expected_copies_have_settled(
            &previous,
            &current,
            diploid(),
            3,
            DEFAULT_CONVERGENCE_THRESHOLD
        ),
        "the reference allele fell by 0.0015 on the frequency scale, half again the 1e-3 \
         threshold, while both alternatives rose by less than it"
    );
}

/// **The loop at three alleles**, because every other fixture here is biallelic and a
/// two-allele locus hides a whole class of defect: the expected copies sum to the
/// cohort's chromosome total, so at two alleles the two movements are equal and opposite
/// and anything that only looks at one of them agrees with the shipped rule by
/// construction (`a_fall_larger_than_every_rise_has_not_settled` is the unit-level
/// counterpart).
///
/// Three samples whose reads pull toward three different homozygotes over three alleles,
/// by the same 2 nats as the biallelic fixture, so the loop takes several passes rather
/// than settling on the first. The asserted property is the one that must hold at any
/// allele count: the cohort's expected copies sum to `ploidy × samples`, because each
/// sample's posterior is normalised over genotypes that each carry `ploidy` copies.
/// Measured, it settles in 2 passes at `[2.4868, 1.7566, 1.7566]`.
#[test]
fn the_loop_settles_at_three_alleles_and_the_copies_still_sum_to_the_chromosomes() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    assert_eq!(view.allele_count(), 3);
    // The table's genotype order at three alleles, read off `genotype_allele_counts`:
    // 0/0, 0/1, 1/1, 0/2, 1/2, 2/2. Each sample favours one homozygote by 2 nats.
    let likelihoods = vec![
        vec![0.0, -2.0, -6.0, -2.0, -6.0, -6.0],
        vec![-6.0, -2.0, 0.0, -6.0, -2.0, -6.0],
        vec![-6.0, -6.0, -6.0, -2.0, -2.0, 0.0],
    ];

    let (outcome, scratch) = loop_over(
        &likelihoods,
        &[1.0, 0.000_5, 0.000_5],
        &alleles,
        &view,
        &RunnableCallingLoopConfig::default(),
    );
    assert!(outcome.converged);

    let copies = scratch.cohort_expected_copies();
    assert_eq!(copies.len(), 3, "one entry per allele");
    let total: f64 = copies.iter().sum();
    assert!(
        (total - 6.0).abs() < 1e-12,
        "three diploid samples carry six chromosomes, and the cohort's expected copies \
         sum to {total} over {copies:?}"
    );
}

/// §6 stops when the movement falls *below* the threshold, and this is the row that sits
/// exactly on it: `2e-3` raw copies over two chromosomes is `1e-3` in `f64` exactly, not
/// a rounding of it, so the strict comparison is what decides the case.
#[test]
fn a_movement_exactly_at_the_threshold_has_not_settled() {
    assert_eq!(
        2e-3_f64 / 2.0,
        DEFAULT_CONVERGENCE_THRESHOLD,
        "the fixture only tests the boundary if it lands on it exactly"
    );
    assert!(
        !cohort_expected_copies_have_settled(
            &[0.0],
            &[2e-3],
            diploid(),
            1,
            DEFAULT_CONVERGENCE_THRESHOLD
        ),
        "§6 stops when the movement falls below the threshold, not when it reaches it"
    );
}

/// A row of no alleles would report **every** locus settled on its first pass, because
/// `all` over nothing is true — the shape an unprepared scratch hands over, and the one
/// door the `NaN` sentinel cannot cover.
#[test]
#[should_panic(expected = "a cohort row of no alleles")]
fn a_convergence_test_over_no_alleles_is_refused() {
    let _ =
        cohort_expected_copies_have_settled(&[], &[], diploid(), 1, DEFAULT_CONVERGENCE_THRESHOLD);
}

/// `zip` stops at the shorter row, so a short previous row would compare the alleles it
/// has and settle on the ones it does not. **A separate assertion from the emptiness
/// check above**, so each test can only be satisfied by the condition it names.
#[test]
#[should_panic(expected = "one entry per allele of the same locus")]
fn cohort_rows_of_different_lengths_are_refused() {
    let _ = cohort_expected_copies_have_settled(
        &[1.0],
        &[1.0, 2.0],
        diploid(),
        1,
        DEFAULT_CONVERGENCE_THRESHOLD,
    );
}

/// A cohort of no samples is zero chromosomes, which turns every movement into an
/// infinity or a `NaN` — so no locus would ever settle and every one would spend the
/// whole cap.
///
/// **The infinite and `NaN` chromosome counts that used to need their own tests are gone
/// rather than untested:** the function takes a [`Ploidy`], which is at least one by
/// construction, and a sample count, so no caller can express them.
#[test]
#[should_panic(expected = "a run whose sample order went missing")]
fn a_convergence_test_over_a_cohort_of_no_samples_is_refused() {
    let _ = cohort_expected_copies_have_settled(
        &[1.0],
        &[1.0],
        diploid(),
        0,
        DEFAULT_CONVERGENCE_THRESHOLD,
    );
}

/// The threshold is a fraction on the frequency scale. `validate` refuses one that is
/// not where a run sets it; this is the check at the point of use, for a threshold built
/// by hand.
#[test]
#[should_panic(expected = "is not one")]
fn a_threshold_that_is_not_a_fraction_is_refused() {
    let _ = cohort_expected_copies_have_settled(&[1.0], &[1.0], diploid(), 1, f64::NAN);
}

/// **The `is_finite()` half of the threshold guard needs an infinity to be reached at
/// all**, because `> 0.0` already refuses zero and `NaN`. And an infinity is the value
/// that matters: it reports every locus settled on its first pass rather than never.
#[test]
#[should_panic(expected = "is not one")]
fn an_infinite_threshold_is_refused() {
    let _ = cohort_expected_copies_have_settled(&[0.0], &[9.0], diploid(), 1, f64::INFINITY);
}

/// **One claimed row per prepared row, in both directions** — the check that replaced the
/// loop's own coefficient-length assertion when the coefficients moved onto the scratch.
///
/// A scratch prepared for three rows and given two samples' coefficients has a row map
/// that describes a different locus: the seeded pass would score two rows and leave the
/// third holding the prior-free initialisation's copies, which the M-step then sums as
/// though this pass had produced them — finite, plausible, wrong and silent.
#[test]
#[should_panic(expected = "of them were claimed")]
fn fewer_claimed_rows_than_the_scratch_was_prepared_for_are_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let (_, _) = loop_over_inbred(
        &three_samples_pulling_apart(),
        &[outbred(), outbred()],
        &[1.0, 0.000_5],
        &alleles,
        &view,
        &RunnableCallingLoopConfig::default(),
    );
}

/// The mirror of the short map, and the same run-shape bug. Without the check it surfaces
/// two modules away as a scratch-indexing panic about the genotype-likelihood table, which
/// names neither the cohort nor the rows.
#[test]
#[should_panic(expected = "of them were claimed")]
fn more_claimed_rows_than_the_scratch_was_prepared_for_are_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let two_samples = three_samples_pulling_apart()[..2].to_vec();
    let (_, _) = loop_over_inbred(
        &two_samples,
        &[outbred(), outbred(), outbred()],
        &[1.0, 0.000_5],
        &alleles,
        &view,
        &RunnableCallingLoopConfig::default(),
    );
}
