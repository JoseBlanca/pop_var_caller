use super::*;

/// **Every intermediate of one sample's E-step, on numbers a reader can follow.**
///
/// A diploid biallelic locus, so three candidate genotypes in the table's order —
/// `0/0` as `[2, 0]`, `0/1` as `[1, 1]`, `1/1` as `[0, 2]`.
///
/// - **the concentration**: seed `[1, 0.5]`, the cohort expected to carry `[3, 1]`
///   copies of which this sample's own are `[1, 1]`, so the other samples showed
///   `[2, 0]` and the concentration is `[3, 0.5]` — exact in binary floating point, so
///   it is asserted exactly.
/// - **the posterior**: likelihoods `ln 1, ln 2, ln 1` and log-priors `ln 1, ln 1,
///   ln 4` give scores whose exponentials are `1, 2, 4`, so the posterior is
///   `1/7, 2/7, 4/7`.
/// - **the expected copies**: the reference allele gets `2·(1/7) + 1·(2/7) = 4/7` and
///   the alternative `1·(2/7) + 2·(4/7) = 10/7`. They sum to 2, the ploidy, because a
///   normalised posterior over genotypes that each carry two copies must.
///
/// **The prior is not decoration here, and the numbers were picked so that dropping it
/// changes the answer loudly.** Measured, by deleting the `+ prior_of_genotype.get()`
/// term: the likelihoods alone give `[0.25, 0.5, 0.25]`, and the test stops at
/// `posterior 0.25 against 0.14285714285714285`.
#[test]
fn one_samples_e_step_matches_the_arithmetic_done_by_hand() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    assert_eq!(view.genotype_count(), 3);
    assert_eq!(view.genotype_allele_counts(), &[2, 0, 1, 1, 0, 2]);

    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[3.0, 1.0]);
    scratch
        .sample_expected_copies_mut(0)
        .copy_from_slice(&[1.0, 1.0]);
    let two = std::f64::consts::LN_2;
    for (slot, value) in scratch
        .sample_genotype_likelihoods_mut(0)
        .iter_mut()
        .zip([0.0, two, 0.0])
    {
        *slot = LogProb(value);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 2.0 * two]),
            inbreeding: outbred(),
        },
    );

    assert_eq!(scratch.sample_concentration(), &[3.0, 0.5]);
    for (got, want) in scratch
        .posterior_row()
        .iter()
        .zip([1.0 / 7.0, 2.0 / 7.0, 4.0 / 7.0])
    {
        assert!((got - want).abs() < 1e-15, "posterior {got} against {want}");
    }
    for (got, want) in scratch
        .sample_expected_copies(0)
        .iter()
        .zip([4.0 / 7.0, 10.0 / 7.0])
    {
        assert!((got - want).abs() < 1e-15, "copies {got} against {want}");
    }
    let total: f64 = scratch.sample_expected_copies(0).iter().sum();
    assert!(
        (total - 2.0).abs() < 1e-15,
        "copies sum to {total}, not the ploidy"
    );
}

/// **The normalisation survives likelihoods no exponential could hold**, because the
/// largest score is subtracted before anything is exponentiated.
///
/// The same shape as the hand-computed case, with every likelihood pushed down by 1,000
/// nats — a locus with about 430 reads at Phred 10 reaches that. The posterior is
/// unchanged, because a constant added to every genotype's score cancels in the
/// normalisation.
///
/// **Measured, with the subtraction removed** (`*slot = slot.exp()`): every term
/// underflows to exactly zero. In debug the normaliser's own check fires first —
/// *"the normaliser cannot come out below one: got 0 over 3 genotypes"* — and under
/// `--release`, where that check is compiled out, the division is `0/0` and the test
/// stops at `posterior NaN against 0.25`. Either way this test fails where the
/// hand-computed one above still passes, which is why it is a test of its own.
#[test]
fn the_posterior_is_unchanged_by_a_constant_that_would_underflow_the_exponentials() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[3.0, 1.0]);
    scratch
        .sample_expected_copies_mut(0)
        .copy_from_slice(&[1.0, 1.0]);
    let two = std::f64::consts::LN_2;
    for (slot, value) in
        scratch
            .sample_genotype_likelihoods_mut(0)
            .iter_mut()
            .zip([-1000.0, -1000.0 + two, -1000.0])
    {
        *slot = LogProb(value);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            inbreeding: outbred(),
        },
    );

    for (got, want) in scratch.posterior_row().iter().zip([0.25, 0.5, 0.25]) {
        assert!((got - want).abs() < 1e-12, "posterior {got} against {want}");
    }
}

/// **At one sample the cohort term is exactly zero, and no line of the E-step tests the
/// cohort size to make it so** (`doc/devel/ng/spec/calling_em_loop.md` §7).
///
/// The cohort's expected copies are that one sample's own, so the leave-one-out
/// subtraction is a number minus itself and the concentration comes back as the seed —
/// asserted with `==` on the floats, because the arithmetic is `seed + (x − x)` and
/// nothing else. This is what makes the loop's fixed point at one sample a consequence
/// rather than a special case.
///
/// **What this test cannot do on its own**, and it is worth knowing which test does:
/// `seed + (x − x) = seed` is a property of [`fill_sample_concentration`], one module
/// away, so any wiring that computes `seed + f(cohort, own)` with `f(x, x) = 0` passes
/// here. Measured: replacing the sample's own copies with the cohort's — which makes
/// every sample's concentration the bare seed at **every** cohort size, so the cohort's
/// evidence never reaches any prior and the loop is inert — leaves this test green, and
/// is killed by `one_samples_e_step_matches_the_arithmetic_done_by_hand`
/// (`left: [1.0, 0.5] right: [3.0, 0.5]`). What this test pins is the *one-sample* claim
/// specifically; the discriminating check on the leave-one-out wiring is that one.
///
/// Run through a shipped prior rather than the stub, because the claim is about what
/// the prior is handed.
#[test]
fn at_one_sample_the_concentration_comes_back_as_the_seed() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    let seed = [1.0, 0.000_5, 0.000_5];
    scratch.seed_concentration_mut().copy_from_slice(&seed);
    let own = [1.7, 0.2, 0.1];
    scratch.cohort_expected_copies_mut().copy_from_slice(&own);
    scratch.sample_expected_copies_mut(0).copy_from_slice(&own);
    for slot in scratch.sample_genotype_likelihoods_mut(0).iter_mut() {
        *slot = LogProb(-3.0);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &MarginalizedDirichletPrior,
            inbreeding: outbred(),
        },
    );

    assert_eq!(scratch.sample_concentration(), &seed);
}

/// **Two properties that hold at every locus, under a shipped prior**: the posterior is
/// a probability distribution over the candidate genotypes, and the expected copies it
/// implies sum to the ploidy.
///
/// **The second is not an independent invariant, and the doc used to say it was.** It
/// follows from the first whenever the fold walks whole genotype rows, because
/// `Σ_a Σ_g p_g·c_{g,a} = Σ_g p_g·(Σ_a c_{g,a})` and every genotype in the table carries
/// exactly the ploidy — so it holds for *any* subset of genotypes whose posterior sums to
/// one, which is why a truncated posterior row passes it while being wrong by a third of
/// a copy. What it does catch is a fold that walks the copy table with the **wrong
/// stride**: measured, folding with `chunks_exact(genotype_count)` instead brings sample
/// 0's copies to **1.9882062059027987** rather than 2 — 0.0118 out, which is why the
/// tolerance here is `1e-12` and not a loose one. The biallelic case above catches that
/// same mutation on a single *value* rather than on a sum: the reference allele still
/// comes out `4/7` and the alternative comes out `0` against `10/7`.
///
/// Every sample here is scored against the same cohort summary but from its own
/// likelihood row, which is also what pins that the function reads the row it was
/// handed: the three samples' posteriors are required to differ.
#[test]
fn every_samples_posterior_is_a_distribution_and_its_copies_sum_to_the_ploidy() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    assert_eq!(view.genotype_count(), 6);
    assert_eq!(view.allele_count(), 3);

    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(3, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.000_5, 0.000_5]);
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[4.4, 1.2, 0.4]);
    let mut first_posterior = Vec::new();
    for sample in 0..3 {
        scratch
            .sample_expected_copies_mut(sample)
            .copy_from_slice(&[1.5, 0.4, 0.1]);
        for (genotype, slot) in scratch
            .sample_genotype_likelihoods_mut(sample)
            .iter_mut()
            .enumerate()
        {
            *slot = LogProb(-(1.0 + sample as f64) * genotype as f64);
        }

        score_one_sample(
            scratch.sample_scoring_buffers_mut(sample),
            &view,
            PassPrior::LeaveOneOut {
                model: &MarginalizedDirichletPrior,
                inbreeding: outbred(),
            },
        );

        let posterior_total: f64 = scratch.posterior_row().iter().sum();
        assert!(
            (posterior_total - 1.0).abs() < 1e-12,
            "sample {sample}'s posterior sums to {posterior_total}"
        );
        assert!(
            scratch
                .posterior_row()
                .iter()
                .all(|p| p.is_finite() && *p >= 0.0),
            "sample {sample}'s posterior is not a distribution: {:?}",
            scratch.posterior_row()
        );
        let copies_total: f64 = scratch.sample_expected_copies(sample).iter().sum();
        assert!(
            (copies_total - 2.0).abs() < 1e-12,
            "sample {sample}'s expected copies sum to {copies_total}, not the ploidy"
        );
        if sample == 0 {
            first_posterior = scratch.posterior_row().to_vec();
        } else {
            assert_ne!(
                first_posterior,
                scratch.posterior_row(),
                "sample {sample} was scored on sample 0's likelihood row"
            );
        }
    }
}

/// **A likelihood row one genotype short is refused, not silently truncated.**
///
/// This is the failure the release-held check exists for: `zip` stops at the shorter
/// side, so without the assertion the last genotype keeps whatever the previous sample's
/// pass left in the posterior buffer and the row is normalised as though it were this
/// sample's. Nothing crashes and the call comes out confident.
///
/// Reached with a hand-built [`SampleScoringBuffers`], because the scratch's own
/// accessor cannot produce a short row — which is the point: the check guards the seam,
/// not the scratch.
#[test]
#[should_panic(expected = "one genotype likelihood per candidate genotype")]
fn a_short_likelihood_row_is_refused() {
    let (_alleles, table) = generic_locus(1);
    let view = table.view();
    let likelihoods = [LogProb(0.0), LogProb(0.0)];
    let mut sample_concentration = vec![0.0; 2];
    let mut workspace = vec![0.0; 2];
    let mut prior_row = vec![LogProb(0.0); 3];
    let mut posterior_row = vec![0.0; 3];
    let mut copies = vec![1.0; 2];

    score_one_sample(
        SampleScoringBuffers {
            sample: 0,
            seed_concentration: &[1.0, 0.5],
            // Consistent with the sample's own `[1.0, 1.0]` below: a cohort total under
            // this sample's own contribution trips a **debug-only** check inside
            // `fill_sample_concentration`, so a fixture with one would keep passing in
            // debug if the length check this test is about were moved or downgraded.
            cohort_expected_copies: &[2.0, 1.0],
            genotype_likelihoods: &likelihoods,
            sample_concentration: &mut sample_concentration,
            prior_per_allele_workspace: &mut workspace,
            prior_row: &mut prior_row,
            posterior_row: &mut posterior_row,
            sample_expected_copies: &mut copies,
        },
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            inbreeding: outbred(),
        },
    );
}

/// **A row with no usable score is raised where it arrives**, rather than turning every
/// posterior entry into a `NaN` that the M-step then sums into the cohort's copies and
/// carries to every other sample's prior on the next pass.
///
/// **This fixture reaches the finiteness check as an `−∞`, not as a `NaN`, and the doc
/// here used to say otherwise.** `largest_score` starts at `−∞` and is only assigned
/// through `score > largest_score`, which every `NaN` loses — so a row that is *entirely*
/// `NaN`, as this one is, never moves the maximum off its starting value and what the
/// check sees is the `−∞`. A `NaN` among finite scores takes the other path, and
/// `a_nan_below_the_largest_score_is_refused` is the test for it.
///
/// Production raised the same condition as an error variant; this module's design makes
/// caller bugs assertions (`doc/devel/ng/spec/calling_em_loop.md` §8).
#[test]
#[should_panic(expected = "genotype scores came out -inf")]
fn a_row_with_no_usable_score_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[2.0, 0.0]);
    scratch
        .sample_expected_copies_mut(0)
        .copy_from_slice(&[2.0, 0.0]);
    for slot in scratch.sample_genotype_likelihoods_mut(0).iter_mut() {
        *slot = LogProb(f64::NAN);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            inbreeding: outbred(),
        },
    );
}
/// **A posterior row one genotype short is refused, not silently truncated.**
///
/// `zip` stops at the shorter side twice over: the score loop writes only as many entries
/// as the row holds, and step 4's fold then walks only that many genotypes. Measured with
/// the check removed, a three-genotype locus scored through a two-entry row gives copies
/// `[1.333…, 0.667…]` against the right `[1.0, 1.0]` — the reference allele over-counted
/// by a third of a copy, with no panic.
///
/// **Neither assertion of the three-sample test notices that**, which is why this test
/// exists rather than leaning on them: the truncated posterior is renormalised to one,
/// and the copies sum to the ploidy for any subset of genotypes whose posterior does.
///
/// Reached with a hand-built [`SampleScoringBuffers`], for the reason
/// `a_short_likelihood_row_is_refused` gives: the scratch cannot produce one.
#[test]
#[should_panic(expected = "one posterior entry per candidate genotype")]
fn a_short_posterior_row_is_refused() {
    let (_alleles, table) = generic_locus(1);
    let view = table.view();
    let likelihoods = [LogProb(0.0), LogProb(std::f64::consts::LN_2), LogProb(0.0)];
    let mut sample_concentration = vec![0.0; 2];
    let mut workspace = vec![0.0; 2];
    let mut prior_row = vec![LogProb(0.0); 3];
    let mut posterior_row = vec![0.0; 2];
    let mut copies = vec![1.0; 2];

    score_one_sample(
        SampleScoringBuffers {
            sample: 0,
            seed_concentration: &[1.0, 0.5],
            cohort_expected_copies: &[3.0, 1.0],
            genotype_likelihoods: &likelihoods,
            sample_concentration: &mut sample_concentration,
            prior_per_allele_workspace: &mut workspace,
            prior_row: &mut prior_row,
            posterior_row: &mut posterior_row,
            sample_expected_copies: &mut copies,
        },
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            inbreeding: outbred(),
        },
    );
}

/// **A seed of the wrong width is refused here, so the panic names the seed.**
///
/// This check is the one of the four that buys a *message* rather than a catch: with it
/// removed, every mis-shaped seed is still refused in release a few lines later, by
/// `fill_sample_concentration` or by `PriorRow::new` — but the message then names one of
/// the buffers that disagrees with the seed rather than the seed itself. The test pins
/// the ordering, so a later step cannot move the check behind those two without noticing.
#[test]
#[should_panic(expected = "one seed concentration per allele")]
fn a_seed_of_the_wrong_width_is_refused() {
    let (_alleles, table) = generic_locus(1);
    let view = table.view();
    let likelihoods = [LogProb(0.0), LogProb(0.0), LogProb(0.0)];
    let mut sample_concentration = vec![0.0; 2];
    let mut workspace = vec![0.0; 2];
    let mut prior_row = vec![LogProb(0.0); 3];
    let mut posterior_row = vec![0.0; 3];
    let mut copies = vec![1.0; 2];

    score_one_sample(
        SampleScoringBuffers {
            sample: 0,
            seed_concentration: &[1.0, 0.5, 0.5],
            cohort_expected_copies: &[2.0, 1.0],
            genotype_likelihoods: &likelihoods,
            sample_concentration: &mut sample_concentration,
            prior_per_allele_workspace: &mut workspace,
            prior_row: &mut prior_row,
            posterior_row: &mut posterior_row,
            sample_expected_copies: &mut copies,
        },
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            inbreeding: outbred(),
        },
    );
}

/// **One `NaN` in a genotype that is not the most probable one is still refused.**
///
/// `largest_score` is only ever assigned through `score > largest_score`, and every
/// comparison against a `NaN` is false, so the maximum is never itself a `NaN` — the
/// finiteness check at the top of the normalisation cannot see this row. Measured under
/// `--release` while the total-weight check was still `debug_assert!`: the call returned
/// normally with every posterior entry and every expected copy `NaN`, which the M-step
/// would then sum into the cohort's copies and carry to every other sample's next prior.
///
/// So the two checks are not redundant, and this is the test that separates them:
/// `a_non_finite_score_is_refused` reaches the finiteness check with an all-`NaN` row
/// that arrives there as the `−∞` the maximum started from, and this one reaches the
/// total-weight check, which is the only `NaN` detector in the function.
#[test]
#[should_panic(expected = "total weight cannot come out below one")]
fn a_nan_below_the_largest_score_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[2.0, 0.0]);
    scratch
        .sample_expected_copies_mut(0)
        .copy_from_slice(&[2.0, 0.0]);
    for (slot, value) in
        scratch
            .sample_genotype_likelihoods_mut(0)
            .iter_mut()
            .zip([f64::NAN, 0.0, 0.0])
    {
        *slot = LogProb(value);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            inbreeding: outbred(),
        },
    );
}

/// **A sample whose own expected copies were never written is refused**, rather than
/// scored against the bare seed with the cohort's evidence silently absent.
///
/// `prepare_for_locus` leaves that row holding `UNWRITTEN_SCRATCH_VALUE`, which is `NaN`,
/// and this is the one place the E-step reads it. **Without a release-held check here the
/// sentinel does not survive to be seen**: measured under `--release` before it was
/// added, the call returned normally and the concentration came back as exactly the seed
/// `[1.0, 0.5]`, because `2.0 − NaN` is `NaN` and `f64::max` returns the *other* operand
/// on a `NaN`, so the leave-one-out term collapses to zero. The two checks that would
/// have caught it — `SampleAlleleCopies::new`'s finiteness check and the count-path
/// desync check — are both `debug_assert!`, one module away.
#[test]
#[should_panic(expected = "a pass reached this sample before anything wrote them")]
fn a_sample_whose_own_copies_were_never_written_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&[2.0, 0.0]);
    for slot in scratch.sample_genotype_likelihoods_mut(0).iter_mut() {
        *slot = LogProb(0.0);
    }
    // `sample_expected_copies` deliberately left as `prepare_for_locus` wrote it.

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            inbreeding: outbred(),
        },
    );
}

/// **At one sample a second pass reproduces the first, bit for bit** — the fixed point
/// `doc/devel/ng/spec/calling_em_loop.md` §13 test 1 names.
///
/// Reachable without the M-step, because at one sample the M-step is `cohort := own`: the
/// cohort's expected copies are that sample's own, so the leave-one-out term stays
/// exactly zero and the prior never moves however many passes run.
///
/// **It cannot stand in for a check on step 4**, and the limit is worth saying rather
/// than leaving for the next reader to find: measured, it still passes with
/// `sample_expected_copies.fill(0.0)` deleted, because whatever step 4 accumulates is
/// copied straight back into the cohort's copies below. It is a check on the
/// concentration path across two passes.
#[test]
fn at_one_sample_a_second_pass_reproduces_the_first() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.000_5, 0.000_5]);
    let own = [1.7, 0.2, 0.1];
    scratch.cohort_expected_copies_mut().copy_from_slice(&own);
    scratch.sample_expected_copies_mut(0).copy_from_slice(&own);
    for (genotype, slot) in scratch
        .sample_genotype_likelihoods_mut(0)
        .iter_mut()
        .enumerate()
    {
        *slot = LogProb(-(genotype as f64));
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &MarginalizedDirichletPrior,
            inbreeding: outbred(),
        },
    );
    let first_posterior = scratch.posterior_row().to_vec();
    let first_copies = scratch.sample_expected_copies(0).to_vec();

    // The M-step at one sample: the cohort's copies are that sample's own.
    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(&first_copies);
    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::LeaveOneOut {
            model: &MarginalizedDirichletPrior,
            inbreeding: outbred(),
        },
    );

    assert_eq!(
        scratch.posterior_row(),
        &first_posterior[..],
        "pass 2's posterior is not pass 1's, so one sample is not a fixed point"
    );
    assert_eq!(
        scratch.sample_expected_copies(0),
        &first_copies[..],
        "pass 2's expected copies are not pass 1's"
    );
}

/// One locus's worth of cohort fixture — what [`run_passes`] scores.
struct CohortFixture<'a> {
    /// Per sample, per candidate genotype: the log-likelihood of that sample's reads.
    likelihoods: &'a [Vec<f64>],
    /// Per sample: the expected allele copies each starts the first pass holding.
    starting_copies: &'a [Vec<f64>],
    /// Which sample of `likelihoods` sits at each index of the run.
    order: &'a [usize],
    /// The cohort row the first pass sees. See [`run_passes`] for why it is given rather
    /// than derived.
    starting_cohort: &'a [f64],
    alleles: &'a CandidateAlleles,
    view: &'a GenotypeTableView<'a>,
    /// The locus's seed concentration.
    seed: &'a [f64],
}

/// Run whole passes over a cohort — the E-step for every sample in index order, then the
/// M-step — with the samples presented in whatever order `order` names.
///
/// Returns each sample's most probable genotype **after the last pass**, in the order it
/// was presented in, and the cohort's summed expected copies.
///
/// **`starting_cohort` is taken rather than derived**, because what it is decides what the
/// counterfactual below even means. Setting it to the *sum* of the samples' own copies
/// gives every sample a leave-one-out term of `n − 1` samples' worth on pass 1; setting it
/// equal to one sample's own copies, where every sample starts alike, makes that term
/// exactly zero — which is the *"seed concentration on its own"* that
/// `spec/calling_em_loop.md` §3 names as the flat pass's only alternative.
///
/// **`flat_first` chooses the first pass's prior**, so the same helper drives the shipped
/// behaviour and the counterfactual `spec/calling_em_loop.md` §3 rules out. A run always
/// passes `true`; `false` exists so a test can show what the other choice converges to.
///
/// **`passes` is the whole reason this helper exists rather than a one-pass one.** On the
/// first pass every sample is scored against a cohort row this function wrote before the
/// loop, so no genotype is downstream of `sum_cohort_expected_copies` and a test reading
/// winners off pass 1 passes against an M-step that computes nothing. Measured in review:
/// with the summing loop deleted, the one-pass version of the permutation
/// test below stayed green while three other tests failed. From pass 2 the winners are
/// scored against the row the previous pass's M-step produced, which is the regime the
/// coupling under test actually occurs in.
///
/// The posterior is one reused buffer, so each sample's winner is read off before the next
/// is scored — which is also why the final pass takes the genotype quality as it scores
/// rather than afterwards.
/// One locus's worth of cohort fixture — what [`run_passes`] scores.
fn run_passes(
    passes: usize,
    flat_first: bool,
    fixture: &CohortFixture<'_>,
) -> (Vec<usize>, Vec<f64>) {
    let &CohortFixture {
        likelihoods,
        starting_copies,
        order,
        starting_cohort,
        alleles,
        view,
        seed,
    } = fixture;
    assert!(passes > 0, "a pass count of zero scores nothing");
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(order.len(), alleles, view);
    scratch.seed_concentration_mut().copy_from_slice(seed);

    scratch
        .cohort_expected_copies_mut()
        .copy_from_slice(starting_cohort);

    for (index, &sample) in order.iter().enumerate() {
        for (slot, &value) in scratch
            .sample_genotype_likelihoods_mut(index)
            .iter_mut()
            .zip(&likelihoods[sample])
        {
            *slot = LogProb(value);
        }
        scratch
            .sample_expected_copies_mut(index)
            .copy_from_slice(&starting_copies[sample]);
    }

    let mut winners = Vec::with_capacity(order.len());
    for pass in 0..passes {
        let prior = if pass == 0 && flat_first {
            PassPrior::Flat
        } else {
            PassPrior::LeaveOneOut {
                model: &MarginalizedDirichletPrior,
                inbreeding: outbred(),
            }
        };
        winners.clear();
        for index in 0..order.len() {
            score_one_sample(scratch.sample_scoring_buffers_mut(index), view, prior);
            let posterior = scratch.posterior_row();
            let (winner, _) = posterior
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.partial_cmp(b).expect("a finite posterior"))
                .expect("a locus has at least one genotype");
            winners.push(winner);
        }
        sum_cohort_expected_copies(scratch.cohort_summing_buffers_mut());
    }
    (winners, scratch.cohort_expected_copies().to_vec())
}

/// The cohort row that is the sum of the samples' own starting copies — what a run's
/// first M-step would have produced from them.
fn cohort_of(starting_copies: &[Vec<f64>], order: &[usize]) -> Vec<f64> {
    let mut cohort = vec![0.0; starting_copies[0].len()];
    for &sample in order {
        for (total, copies) in cohort.iter_mut().zip(&starting_copies[sample]) {
            *total += copies;
        }
    }
    cohort
}

/// Three samples whose reads pull them to three different genotypes, with different
/// starting copies — the fixture the two cohort-wide tests share.
fn three_disagreeing_samples() -> (Vec<Vec<f64>>, Vec<Vec<f64>>, [f64; 3]) {
    let likelihoods = vec![
        vec![0.0, -9.0, -18.0, -9.0, -18.0, -18.0],
        vec![-9.0, 0.0, -9.0, -9.0, -9.0, -18.0],
        vec![-18.0, -9.0, 0.0, -9.0, -9.0, -18.0],
    ];
    let starting_copies = vec![
        vec![1.9, 0.05, 0.05],
        vec![0.6, 1.3, 0.1],
        vec![1.2, 0.2, 0.6],
    ];
    (likelihoods, starting_copies, [1.0, 0.000_5, 0.000_5])
}

/// **Presenting the cohort's samples in a different order calls the same genotypes.**
///
/// The M-step's inputs are a set, not a sequence: which sample sits at which index is an
/// accident of how the run was assembled, and no call may depend on it. Three samples with
/// different reads and different starting copies are run in `[0, 1, 2]` and again in
/// `[2, 0, 1]`, and each sample's most probable genotype is required to follow it.
///
/// **Two passes, and the second is the only one that tests anything here.** On pass 1
/// every sample is scored against a cohort row the fixture wrote, so no genotype depends
/// on the M-step at all; the winners are taken after pass 2, when they are scored against
/// the row pass 1's M-step produced. Measured: written as one pass, this test passed with
/// the summing loop **deleted entirely**, while three other tests failed.
///
/// **All three samples must call different genotypes**, and the assertion below insists on
/// it rather than on the weaker "some two differ": if two samples agree, an implementation
/// that swapped exactly those two satisfies every assertion here. The fixture's earlier
/// version called `[0, 1, 0]` and passed a guard that allowed it.
///
/// **The cohort's summed copies are compared with a tolerance and the genotypes are not**,
/// and that difference is the whole of `spec/calling_em_loop.md` §13 test 2. Reordering a
/// floating-point sum moves it in the last bits, so an exact comparison here would fail on
/// a correct implementation; the genotypes must not move at all. The *bitwise* check
/// belongs to a fixture built for it — [`the_sum_runs_in_ascending_sample_order`].
#[test]
fn presenting_the_samples_in_another_order_calls_the_same_genotypes() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    let (likelihoods, starting_copies, seed) = three_disagreeing_samples();

    let (called_in_order, cohort_in_order) = run_passes(
        2,
        true,
        &CohortFixture {
            likelihoods: &likelihoods,
            starting_copies: &starting_copies,
            order: &[0, 1, 2],
            starting_cohort: &cohort_of(&starting_copies, &[0, 1, 2]),
            alleles: &alleles,
            view: &view,
            seed: &seed,
        },
    );
    let (called_rotated, cohort_rotated) = run_passes(
        2,
        true,
        &CohortFixture {
            likelihoods: &likelihoods,
            starting_copies: &starting_copies,
            order: &[2, 0, 1],
            starting_cohort: &cohort_of(&starting_copies, &[2, 0, 1]),
            alleles: &alleles,
            view: &view,
            seed: &seed,
        },
    );

    // The fixture is only a test of the ordering if **all three** samples disagree; two
    // that agree are two this test cannot see swapped.
    assert_eq!(called_in_order.len(), 3);
    assert!(
        called_in_order[0] != called_in_order[1]
            && called_in_order[1] != called_in_order[2]
            && called_in_order[0] != called_in_order[2],
        "two samples called the same genotype, so this fixture cannot see them swapped: \
         {called_in_order:?}"
    );
    for (presented, &sample) in [2usize, 0, 1].iter().enumerate() {
        assert_eq!(
            called_rotated[presented], called_in_order[sample],
            "sample {sample} was called {} presented in order and {} presented at index \
             {presented}",
            called_in_order[sample], called_rotated[presented]
        );
    }
    for (allele, (rotated_total, in_order_total)) in
        cohort_rotated.iter().zip(&cohort_in_order).enumerate()
    {
        assert!(
            (rotated_total - in_order_total).abs() < 1e-12,
            "allele {allele}'s cohort copies moved from {in_order_total} to \
             {rotated_total} on a reordering"
        );
    }
}

/// **The sum runs in ascending sample order, and the check is on the bits.**
///
/// This is the mutation oracle `spec/calling_em_loop.md` §13 test 2 asks for. The reference
/// allele's column is four samples carrying `2⁻⁵³`, `1.0`, `2⁻⁵³` and `2⁻⁵²` copies —
/// `2⁻⁵³` being exactly half the gap between `1.0` and its neighbour, so it rounds away
/// against the one and does not round away against another `2⁻⁵³`.
///
/// **The fixture certifies itself.** Rather than quoting what the wrong orders give, the
/// test builds them and asserts they differ: reversal, and every adjacent transposition
/// **except the first**. That exception is not a gap — `t = a; t += b` and `t = b; t += a`
/// are bit-identical for every pair, IEEE addition being commutative, so **no fixture of
/// any shape can separate a swap of the first two samples, and no implementation can be
/// wrong by making one.** An earlier three-sample fixture separated reversal only, and a
/// walk that swapped rows 1 and 2 passed it.
///
/// **What it still cannot see**: 3 of the 23 non-identity permutations of four samples sum
/// bit-identically to ascending. Floating-point addition is commutative in pairs, so no
/// column separates every permutation; this one separates the ones an implementation could
/// plausibly produce.
///
/// The cost of getting this wrong is not a wrong answer at one locus; it is a run whose
/// output depends on the worker count (`spec/calling_em_loop.md` §8).
#[test]
fn the_sum_runs_in_ascending_sample_order() {
    let half_an_ulp_at_one = f64::from_bits(0x3CA0_0000_0000_0000); // 2^-53
    let one_ulp_at_one = f64::from_bits(0x3CB0_0000_0000_0000); // 2^-52
    assert_eq!(half_an_ulp_at_one, float::powi(2.0, -53));
    assert_eq!(one_ulp_at_one, float::powi(2.0, -52));

    // Sample-major, two alleles: the reference column is what carries the property, and
    // the alternative column is three-and-a-bit copies of an order-invariant `1.0`.
    let reference = [half_an_ulp_at_one, 1.0, half_an_ulp_at_one, one_ulp_at_one];
    let table: Vec<f64> = reference.iter().flat_map(|&r| [r, 1.0]).collect();

    let fold_in = |walk: &[usize]| -> f64 {
        walk.iter()
            .skip(1)
            .fold(reference[walk[0]], |total, &row| total + reference[row])
    };
    let ascending = fold_in(&[0, 1, 2, 3]);

    // Reversal, and every adjacent transposition but the first, must be visible. The first
    // cannot be: `a + b` and `b + a` are the same bits.
    assert_ne!(
        fold_in(&[3, 2, 1, 0]).to_bits(),
        ascending.to_bits(),
        "this fixture cannot see a reversed sum"
    );
    for (i, walk) in [[1, 0, 2, 3], [0, 2, 1, 3], [0, 1, 3, 2]]
        .iter()
        .enumerate()
    {
        let swapped = fold_in(walk);
        if i == 0 {
            assert_eq!(
                swapped.to_bits(),
                ascending.to_bits(),
                "swapping the first two samples must be a no-op, IEEE addition being \
                 commutative"
            );
        } else {
            assert_ne!(
                swapped.to_bits(),
                ascending.to_bits(),
                "this fixture cannot see samples {i} and {} swapped",
                i + 1
            );
        }
    }

    let mut cohort = vec![f64::NAN; 2];
    sum_cohort_expected_copies(CohortSummingBuffers {
        sample_count: 4,
        per_sample_expected_copies: &table,
        cohort_expected_copies: &mut cohort,
    });

    assert_eq!(
        cohort[0].to_bits(),
        ascending.to_bits(),
        "the reference allele's copies came to {} where ascending sample order gives {}",
        cohort[0],
        ascending
    );
    assert_eq!(cohort[1], 4.0);
}

/// **The cohort's copies come to the ploidy times the sample count**, because every
/// sample's own copies come to the ploidy and the M-step only adds them up.
///
/// It is the cheapest check there is on a sum that dropped or double-counted a sample —
/// either moves the total by a whole 2.0, so the `1e-12` tolerance is ample.
///
/// **The grand total is asserted and so is each allele**, because the total alone cannot
/// see an M-step that mixed the alleles up: any permutation of the cohort row sums to the
/// same number. The per-allele half compares against a fold written independently here.
///
/// **The rows are the samples the loop scores.** A sample the candidate step set aside
/// (`spec/calling_em_loop.md` §5.0) gets no row, so at a real locus the total is the ploidy
/// times the *callable* count.
#[test]
fn the_cohort_carries_the_ploidy_for_every_sample_it_summed() {
    let (alleles, table) = generic_locus(2);
    let view = table.view();
    let (likelihoods, starting_copies, seed) = three_disagreeing_samples();

    let (_winners, cohort) = run_passes(
        1,
        false,
        &CohortFixture {
            likelihoods: &likelihoods,
            starting_copies: &starting_copies,
            order: &[0, 1, 2],
            starting_cohort: &cohort_of(&starting_copies, &[0, 1, 2]),
            alleles: &alleles,
            view: &view,
            seed: &seed,
        },
    );

    let total: f64 = cohort.iter().sum();
    let expected = f64::from(view.ploidy().get()) * 3.0;
    assert!(
        (total - expected).abs() < 1e-12,
        "the cohort carries {total} copies over three diploid samples, not {expected}"
    );

    // Per allele too, because the total alone cannot see a cohort row whose alleles were
    // permuted: any permutation sums to the same number.
    assert_eq!(cohort.len(), view.allele_count());
    assert!(
        cohort.iter().all(|c| c.is_finite() && *c >= 0.0),
        "every cohort entry is a count of genome copies: {cohort:?}"
    );
    // The reference allele must carry more than either alternative here: all three
    // samples start reference-heavy and only one is pulled off it hard.
    assert!(
        cohort[0] > cohort[1] && cohort[0] > cohort[2],
        "the reference allele should dominate this fixture's cohort: {cohort:?}"
    );
}

/// **A sample row the E-step never wrote is refused, and one check over the alleles is
/// enough to catch it** — `NaN` survives every addition, so an unwritten row anywhere in a
/// table of any size reaches the cohort entry it belongs to.
///
/// This is what stands between the scratch's sentinel and a cohort summary that every
/// sample's next prior is built from.
#[test]
#[should_panic(expected = "a sample row the E-step never wrote")]
fn a_sample_row_that_was_never_written_is_refused() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &view);
    scratch
        .sample_expected_copies_mut(0)
        .copy_from_slice(&[1.0, 1.0]);
    // Sample 1's row deliberately left as `prepare_for_locus` wrote it.

    sum_cohort_expected_copies(scratch.cohort_summing_buffers_mut());
}

/// **A per-sample table that is not `samples × alleles` is refused**, rather than summed
/// over whatever rows `chunks_exact` happens to find. A table one row short would sum a
/// cohort of `n − 1` and report it as `n`, which no later check can see: the copies would
/// be finite, non-negative, and smaller than they should be by one sample's worth.
#[test]
#[should_panic(expected = "the per-sample copies are samples × alleles")]
fn a_per_sample_table_of_the_wrong_size_is_refused() {
    let mut cohort = vec![f64::NAN; 2];
    sum_cohort_expected_copies(CohortSummingBuffers {
        sample_count: 3,
        per_sample_expected_copies: &[1.0, 1.0, 1.0, 1.0],
        cohort_expected_copies: &mut cohort,
    });
}
/// **A cohort of no samples is refused, and the size check cannot do it.**
///
/// A count of zero *satisfies* the shape check, because `0 == 0 × alleles` — so without
/// this guard the walk sums nothing and hands back a cohort row that reads like a locus
/// where nobody carries anything, which the finiteness check accepts. Measured in review:
/// with both `> 0` guards deleted the whole of `calling` stays green.
#[test]
#[should_panic(expected = "a sum over none is a run whose sample order went missing")]
fn a_cohort_of_no_samples_is_refused() {
    let mut cohort = vec![f64::NAN; 2];
    sum_cohort_expected_copies(CohortSummingBuffers {
        sample_count: 0,
        per_sample_expected_copies: &[],
        cohort_expected_copies: &mut cohort,
    });
}

/// **A locus of no alleles is refused**, for the same reason and by the same reasoning:
/// `entries == samples × 0` holds for any sample count, so the shape check passes and the
/// walk writes nothing.
#[test]
#[should_panic(expected = "a cohort row of no alleles")]
fn a_locus_of_no_alleles_is_refused() {
    let mut cohort: Vec<f64> = Vec::new();
    sum_cohort_expected_copies(CohortSummingBuffers {
        sample_count: 3,
        per_sample_expected_copies: &[],
        cohort_expected_copies: &mut cohort,
    });
}

/// **An infinite total is refused** — the half of the output check that `>= 0.0` cannot
/// do, and the one a real overflow would reach.
///
/// The two halves are not redundant and each needs its own case:
/// [`a_negative_total_is_refused`] is the other.
#[test]
#[should_panic(expected = "must be finite and at or above zero")]
fn an_infinite_total_is_refused() {
    let mut cohort = vec![f64::NAN; 2];
    sum_cohort_expected_copies(CohortSummingBuffers {
        sample_count: 2,
        per_sample_expected_copies: &[f64::MAX, 1.0, f64::MAX, 1.0],
        cohort_expected_copies: &mut cohort,
    });
}

/// **A negative total is refused** — the half of the output check that `is_finite` cannot
/// do.
///
/// **And this test is also where the check's limit lives.** It catches a total that comes
/// out negative; it does **not** catch negative *inputs* that cancel, because two rows of
/// `-1.0` and `3.0` sum to a perfectly acceptable `2.0`. Checking every input instead
/// would cost `samples × alleles` on the loop's hot path to catch a state no producer can
/// reach, expected copies being a posterior-weighted sum of non-negative counts. Reached
/// with a hand-built bundle for that reason: the E-step cannot produce one.
#[test]
#[should_panic(expected = "must be finite and at or above zero")]
fn a_negative_total_is_refused() {
    let mut cohort = vec![f64::NAN; 2];
    sum_cohort_expected_copies(CohortSummingBuffers {
        sample_count: 2,
        per_sample_expected_copies: &[-3.0, 1.0, 1.0, 1.0],
        cohort_expected_copies: &mut cohort,
    });
}

/// **At one sample the cohort's copies are that sample's own row, bit for bit.**
///
/// The boundary the whole design leans on: `spec/calling_em_loop.md` §7 makes one sample a
/// first-class case rather than a degraded one, and the M-step there is not an
/// approximation of a sum — it *is* the row. Asserted on the bits, because a sum of one
/// addend has no rounding to hide behind, so anything that scales, doubles or reorders the
/// addends shows up.
///
/// Nothing else in this module's tests runs the M-step at one sample.
#[test]
fn the_cohort_of_one_sample_is_that_samples_own_row_bit_for_bit() {
    let own = [0.7_f64, 1.3];
    let mut cohort = vec![f64::NAN; 2];
    sum_cohort_expected_copies(CohortSummingBuffers {
        sample_count: 1,
        per_sample_expected_copies: &own,
        cohort_expected_copies: &mut cohort,
    });
    assert_eq!(cohort[0].to_bits(), own[0].to_bits());
    assert_eq!(cohort[1].to_bits(), own[1].to_bits());
}

/// **From pass 2 the `NaN` sentinel is gone, and a row this pass did not write is summed
/// as though it were current.** This is the limit of the finiteness check, pinned so that
/// the doc comment claiming it and the behaviour cannot drift apart.
///
/// The scratch fills the per-sample rows with `UNWRITTEN_SCRATCH_VALUE` once per locus,
/// not once per pass — deliberately, because `score_one_sample` reads each sample's
/// previous copies as its leave-one-out term. So a sample whose E-step is skipped on a
/// later pass contributes its **previous** copies, finite and plausible, and nothing
/// raises.
///
/// **This is not a bug in the M-step; it is why a loop that skipped a sample mid-run would
/// need a per-pass written mask.** The driver never does — a sample set aside under spec
/// §5.0 has no row from the start — so this pins a hazard for any future caller rather than
/// a failure a run can reach.
#[test]
fn from_the_second_pass_an_unwritten_row_is_summed_as_current() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(2, &alleles, &view);

    // Pass 1: both samples write their rows, and the sum is what they say.
    scratch
        .sample_expected_copies_mut(0)
        .copy_from_slice(&[2.0, 0.0]);
    scratch
        .sample_expected_copies_mut(1)
        .copy_from_slice(&[0.0, 2.0]);
    sum_cohort_expected_copies(scratch.cohort_summing_buffers_mut());
    assert_eq!(scratch.cohort_expected_copies(), &[2.0, 2.0]);

    // Pass 2: only sample 0 writes. Sample 1's pass-1 row is still there.
    scratch
        .sample_expected_copies_mut(0)
        .copy_from_slice(&[1.0, 1.0]);
    sum_cohort_expected_copies(scratch.cohort_summing_buffers_mut());
    assert_eq!(
        scratch.cohort_expected_copies(),
        &[1.0, 3.0],
        "sample 1's pass-1 copies were summed into pass 2's cohort, which is the limit \
         the finiteness check cannot see"
    );
}

proptest! {
    /// **The sum is an ascending fold, at every shape and every value** — the property the
    /// five hand-built fixtures can only sample.
    ///
    /// This function exists *because* floating-point addition is not associative, so its
    /// contract is not "the total" but one specific fold in one specific order. The
    /// expectation here is written independently of the implementation, as a plain
    /// ascending fold over the same table, and compared **on the bits**.
    #[test]
    fn the_sum_is_an_ascending_fold_at_every_shape(
        sample_count in 1_usize..24,
        allele_count in 1_usize..8,
        raw in prop::collection::vec(0.0_f64..2.0, 24 * 8),
    ) {
        let table: Vec<f64> = raw[..sample_count * allele_count].to_vec();
        let mut cohort = vec![f64::NAN; allele_count];
        sum_cohort_expected_copies(CohortSummingBuffers {
            sample_count,
            per_sample_expected_copies: &table,
            cohort_expected_copies: &mut cohort,
        });

        for allele in 0..allele_count {
            let expected = (1..sample_count).fold(table[allele], |total, sample| {
                total + table[sample * allele_count + allele]
            });
            prop_assert_eq!(
                cohort[allele].to_bits(),
                expected.to_bits(),
                "allele {} over {} samples", allele, sample_count
            );
        }
    }
}
/// **The flat first pass runs at a locus whose cohort summary does not exist yet**, which
/// is the whole reason it is a variant of this function rather than a flat prior handed to
/// it.
///
/// Neither the cohort's expected copies nor the sample's own are written here — they hold
/// the scratch's `NaN` sentinel, exactly as `prepare_for_locus` left them, because on the
/// first pass through a locus no pass has produced them. A seeded pass would read them:
/// **measured, handing this same scratch `PassPrior::LeaveOneOut` panics in debug** with
/// *"the cohort's expected allele copies … must be finite"*, and in release would return a
/// concentration equal to the bare seed, since `f64::max` returns the other operand on a
/// `NaN`. That is the seeded first pass `spec/calling_em_loop.md` §3 rules out.
///
/// The posterior must be the likelihoods alone: `ln 1, ln 2, ln 1` normalise to
/// `0.25, 0.5, 0.25`, and the copies to one of each allele.
#[test]
fn the_flat_first_pass_needs_no_cohort_summary() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    // The seed is filled because it is a property of the locus; the cohort's copies and
    // the sample's own deliberately are not.
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    assert!(
        scratch.cohort_expected_copies().iter().all(|c| c.is_nan()),
        "the fixture is only a test of the flat pass if the cohort summary is unwritten"
    );
    let two = std::f64::consts::LN_2;
    for (slot, value) in scratch
        .sample_genotype_likelihoods_mut(0)
        .iter_mut()
        .zip([0.0, two, 0.0])
    {
        *slot = LogProb(value);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::Flat,
    );

    for (got, want) in scratch.posterior_row().iter().zip([0.25, 0.5, 0.25]) {
        assert!((got - want).abs() < 1e-15, "posterior {got} against {want}");
    }
    for (got, want) in scratch.sample_expected_copies(0).iter().zip([1.0, 1.0]) {
        assert!((got - want).abs() < 1e-15, "copies {got} against {want}");
    }
}

/// **A flat first pass finds the variant on pass 1; a seeded one takes nine passes to get
/// there.** This is `spec/calling_em_loop.md` §13 test 3's trap — and the measurement is
/// weaker than §3's account of it, which is why the numbers are here rather than the
/// story.
///
/// A thin locus: every sample's reads favour the heterozygote over the homozygous
/// reference by **1 nat**, about 4.3 Phred — one alternative read among a handful. The seed
/// is `[1.0, 0.000_5]`, under which the homozygous reference outweighs the heterozygote by
/// **7.6009 nats** (measured: the row prints `[0.693, −6.908, −7.600]` at that
/// concentration), so on a seeded first pass the reads lose by 6.6 nats.
///
/// **What is measured, at 63 samples**, alternative copies per sample:
///
/// | pass | flat start | seeded start |
/// |---|---|---|
/// | 1 | 0.731 — every sample heterozygous | 0.0014 — every sample homozygous reference |
/// | 6 | 0.767 | 0.151 — still homozygous reference |
/// | 9 | 0.767 | 0.633 — **flips to heterozygous** |
/// | 30 | 0.767332 | 0.767332 |
///
/// **So both starts reach the same answer, and what the flat pass buys on this fixture is
/// eight passes rather than a different call.** `spec/calling_em_loop.md` §3 says the
/// seeded loop *"converges, and it converges to no-variant, having never let the reads
/// speak"* — that stronger claim is **not** what this fixture shows, and no fixture built
/// for this test showed it: a rare-variant shape, where a handful of carriers sit among 60
/// samples whose reads are firmly homozygous reference, was swept over carrier count and
/// likelihood advantage and the two starts agreed in **every** cell. Whether the delay ever
/// becomes permanent is governed by the pass cap, which ships at 50, and by
/// §12's question 7 on real data.
///
/// **The delay is not nothing.** Production's own comment records the expectation-maximization
/// converging in 3 to 5 passes, so a locus that needs 9 under one start and 1 under the
/// other is a locus whose answer depends on where the cap falls.
///
/// **Where the two starts part depends on both the reads and the cohort size**, swept over
/// likelihood advantage × samples (3, 6, 20, 63): at **2 nats and above they agree at every
/// size** — the reads beat the seed either way and the choice costs nothing. At **0.5 nats**
/// they part at 20 and 63 samples but not at 3 or 6, where the flat start does not reach the
/// heterozygote either. At **1 nat** they part at every size tried, which is why the tests
/// use it. **An earlier version of this comment reported that sweep as though cohort size
/// were not an axis**, which is the failure `CLAUDE.md` names — a figure measured in one
/// corner and reported as a property of the caller.
///
/// Run at 3 samples and at 63 — the tomato panel's size — because a difference visible at
/// one cohort size would be an artifact of that size.
#[test]
fn a_flat_first_pass_finds_the_variant_at_once_where_a_seeded_one_takes_nine_passes() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    assert_eq!(view.genotype_count(), 3, "0/0, 0/1, 1/1");
    let seed = [1.0, 0.000_5];

    for samples in [3usize, 63] {
        let likelihoods = vec![vec![-1.0, 0.0, -30.0]; samples];
        let starting_copies = vec![vec![1.0, 1.0]; samples];
        let order: Vec<usize> = (0..samples).collect();
        let fixture = CohortFixture {
            likelihoods: &likelihoods,
            starting_copies: &starting_copies,
            order: &order,
            starting_cohort: &starting_copies[0],
            alleles: &alleles,
            view: &view,
            seed: &seed,
        };

        // One pass is all the flat start needs.
        let (flat_first_pass, flat_cohort) = run_passes(1, true, &fixture);
        assert!(
            flat_first_pass.iter().all(|&g| g == 1),
            "at {samples} samples the flat start should call every sample heterozygous \
             on pass 1 (genotype 1 of 0/0, 0/1, 1/1): {flat_first_pass:?}"
        );
        assert!(
            flat_cohort[1] / (samples as f64) > 0.7,
            "the flat start's first pass should leave about 0.73 copies of the \
             alternative per sample: {flat_cohort:?}"
        );

        // Six passes and the seeded start is still calling no-variant.
        let (seeded_at_six, seeded_cohort) = run_passes(6, false, &fixture);
        assert!(
            seeded_at_six.iter().all(|&g| g == 0),
            "at {samples} samples the seeded start should still call every sample \
             homozygous reference at six passes: {seeded_at_six:?}"
        );
        assert!(
            seeded_cohort[1] / (samples as f64) < 0.2,
            "at six passes the seeded start should still have the alternative all but \
             absent: {seeded_cohort:?}"
        );

        // And by thirty they agree — the delay, not a different answer.
        let (_flat_late, flat_late_cohort) = run_passes(30, true, &fixture);
        let (_seeded_late, seeded_late_cohort) = run_passes(30, false, &fixture);
        assert!(
            (flat_late_cohort[1] - seeded_late_cohort[1]).abs() < 1e-4,
            "at thirty passes the two starts should have reached the same cohort \
             frequency at {samples} samples: {flat_late_cohort:?} against \
             {seeded_late_cohort:?}"
        );
    }
}

/// **After the flat pass the expected copies reflect the reads and not the seed** —
/// `spec/calling_em_loop.md` §13 test 3's first half, on the same fixture as the trap above
/// so the two read together.
///
/// The seed puts about 2,000 to 1 on the homozygous reference. One flat pass over reads
/// that favour the heterozygote by a single nat leaves the cohort carrying more than half a
/// copy of the alternative allele per sample — which is the number the second pass's prior
/// is built from, and the whole point of starting flat.
#[test]
fn after_the_flat_pass_the_copies_reflect_the_reads_and_not_the_seed() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let samples = 6;
    let likelihoods = vec![vec![-1.0, 0.0, -30.0]; samples];
    let starting_copies = vec![vec![1.0, 1.0]; samples];
    let order: Vec<usize> = (0..samples).collect();

    let (_winners, cohort) = run_passes(
        1,
        true,
        &CohortFixture {
            likelihoods: &likelihoods,
            starting_copies: &starting_copies,
            order: &order,
            starting_cohort: &starting_copies[0],
            alleles: &alleles,
            view: &view,
            seed: &[1.0, 0.000_5],
        },
    );

    // Six diploid samples carry twelve copies in all, and the reads put more than half of
    // them on the alternative allele — where the seed alone would leave near none.
    let total: f64 = cohort.iter().sum();
    assert!(
        (total - 12.0).abs() < 1e-12,
        "twelve copies over six diploid samples: {total}"
    );
    assert!(
        cohort[1] / samples as f64 > 0.5,
        "the reads should put more than half a copy of the alternative on each sample, \
         not the seed's near-zero: {cohort:?}"
    );
}
/// **The flat pass leaves the cohort's copies untouched, and this is the only test that
/// says so in the release profile.**
///
/// It is the property the whole variant exists for, and until this test the only thing
/// standing behind it was a `debug_assert` inside `SampleAlleleCopies::new` one module
/// away — so a mutant that built the concentration anyway passed the release CI step.
/// Asserted positively, on the buffers, rather than by not panicking: after a flat pass the
/// cohort row and the sample's concentration are both still the `NaN` sentinel
/// `prepare_for_locus` wrote, because nothing read or wrote either.
///
/// Measured: under a mutant that builds the concentration on the flat arm, the
/// concentration comes back `[1.0, 0.5]` — the bare seed — where it is `[NaN, NaN]` as
/// shipped.
#[test]
fn the_flat_pass_touches_neither_the_cohort_summary_nor_the_concentration() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.5]);
    for slot in scratch.sample_genotype_likelihoods_mut(0).iter_mut() {
        *slot = LogProb(0.0);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::Flat,
    );

    assert!(
        scratch.cohort_expected_copies().iter().all(|c| c.is_nan()),
        "a flat pass must not write the cohort's copies: {:?}",
        scratch.cohort_expected_copies()
    );
    assert!(
        scratch.sample_concentration().iter().all(|c| c.is_nan()),
        "a flat pass builds no concentration, so the buffer must still hold the sentinel: \
         {:?}",
        scratch.sample_concentration()
    );
}

/// **A prior row of the wrong width is refused on the flat arm too.**
///
/// This width was once checked only inside `PriorRow::new`, which a flat pass never
/// calls — so the seeded arm panicked in release and the flat arm returned a **wrong
/// posterior in silence**: measured at a 3-genotype locus given a 2-entry prior row,
/// `[0.199, 0.399, 0.402]` against the right `[0.25, 0.5, 0.25]`, the tail entry being the
/// stale value the buffer arrived holding.
#[test]
#[should_panic(expected = "one prior entry per candidate genotype")]
fn a_prior_row_of_the_wrong_width_is_refused_on_a_flat_pass() {
    let (_alleles, table) = generic_locus(1);
    let view = table.view();
    let likelihoods = [LogProb(0.0), LogProb(0.0), LogProb(0.0)];
    let mut sample_concentration = vec![0.0; 2];
    let mut workspace = vec![0.0; 2];
    let mut prior_row = vec![LogProb(0.0); 2]; // one short
    let mut posterior_row = vec![0.0; 3];
    let mut copies = vec![1.0; 2];

    score_one_sample(
        SampleScoringBuffers {
            sample: 0,
            seed_concentration: &[1.0, 0.5],
            cohort_expected_copies: &[2.0, 1.0],
            genotype_likelihoods: &likelihoods,
            sample_concentration: &mut sample_concentration,
            prior_per_allele_workspace: &mut workspace,
            prior_row: &mut prior_row,
            posterior_row: &mut posterior_row,
            sample_expected_copies: &mut copies,
        },
        &view,
        PassPrior::Flat,
    );
}

/// **An expected-copies row of the wrong width is refused on the flat arm too.**
///
/// Its width was `fill_sample_concentration`'s to check, and a flat pass never calls it.
/// Measured at a 3-allele locus given a 2-entry copies row, the flat arm returned
/// `[0.667, 0.667]` — a diploid sample carrying 1.33 copies of a genome — with nothing
/// raised.
#[test]
#[should_panic(expected = "one expected-copies entry per allele")]
fn an_expected_copies_row_of_the_wrong_width_is_refused_on_a_flat_pass() {
    let (_alleles, table) = generic_locus(2);
    let view = table.view();
    let likelihoods = vec![LogProb(0.0); 6];
    let mut sample_concentration = vec![0.0; 3];
    let mut workspace = vec![0.0; 3];
    let mut prior_row = vec![LogProb(0.0); 6];
    let mut posterior_row = vec![0.0; 6];
    let mut copies = vec![1.0; 2]; // one short

    score_one_sample(
        SampleScoringBuffers {
            sample: 0,
            seed_concentration: &[1.0, 0.5, 0.5],
            cohort_expected_copies: &[2.0, 1.0, 1.0],
            genotype_likelihoods: &likelihoods,
            sample_concentration: &mut sample_concentration,
            prior_per_allele_workspace: &mut workspace,
            prior_row: &mut prior_row,
            posterior_row: &mut posterior_row,
            sample_expected_copies: &mut copies,
        },
        &view,
        PassPrior::Flat,
    );
}

/// **A sample with no reads puts no copies in the cohort's sum, even on the flat pass**
/// (`spec/calling_em_loop.md` §5.0.2).
///
/// On a flat pass there is no prior, so its posterior is the normalised genotype table and
/// its copies would be the average genotype: a full copy of the alternative allele at a
/// biallelic locus, the same vote as a confident heterozygote. This test used to pin that
/// `1.0`, raised for the owner as undecided; the ruling of 2026-10-31 is that a sample
/// whose reads said nothing contributes nothing, on every pass.
#[test]
fn a_sample_with_no_reads_contributes_no_copies_to_the_flat_pass() {
    let (alleles, table) = generic_locus(1);
    let view = table.view();
    let mut scratch = CallingScratch::<()>::default();
    scratch.prepare_for_locus(1, &alleles, &view);
    scratch
        .seed_concentration_mut()
        .copy_from_slice(&[1.0, 0.000_5]);
    // A sample with no reads: every genotype equally likely.
    for slot in scratch.sample_genotype_likelihoods_mut(0).iter_mut() {
        *slot = LogProb(0.0);
    }

    score_one_sample(
        scratch.sample_scoring_buffers_mut(0),
        &view,
        PassPrior::Flat,
    );

    assert_eq!(scratch.sample_expected_copies(0), &[0.0, 0.0]);
}
