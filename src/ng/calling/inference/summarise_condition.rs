//! **Summarise the cohort, then condition each sample on the summary** — the default of the
//! three ways this caller can score a cohort, and the one the rest of the design describes.
//!
//! Rather than score every sample's genotype jointly with every other sample's, it replaces
//! the other samples with a summary of them — how many copies of each allele the cohort is
//! expected to carry at this locus — and scores each sample against that
//! (`doc/devel/ng/spec/calling_em_loop.md` §2). The architecture calls this *arm A*; the
//! other two score whole cohort-wide genotype assignments jointly and are not built
//! (`doc/devel/ng/arch/calling_em_loop.md` §3.1).
//!
//! **The summary and the scores are each other's inputs**, which is why this is a loop: the
//! copies come from the scores and the scores come from the copies. One turn of it has two
//! halves, and the design names them after the algorithm it inherits,
//! expectation-maximization:
//!
//! - the **E-step** — for each sample, how probable each candidate genotype is, given that
//!   sample's reads and the summary of the others. [`score_one_sample`] is one sample's
//!   share of it, and this step is where all of the loop's arithmetic lives.
//! - the **M-step** — the summary again, from the scores the E-step just produced: the
//!   cohort's expected allele copies, summed over the samples.
//!
//! **Nothing here allocates.** Every buffer belongs to the worker's
//! [`CallingScratch`](crate::ng::calling::CallingScratch) and is reused at every locus and
//! on every pass, which is why the arithmetic below takes a
//! [`SampleScoringBuffers`] rather than returning anything.

use crate::ng::calling::GenotypeTableView;
use crate::ng::calling::SampleScoringBuffers;
use crate::ng::calling::genotype_prior::{
    CohortAlleleCopies, Concentration, GenotypePriorModel, PriorRow, SampleAlleleCopies,
    fill_sample_concentration,
};
use crate::ng::types::InbreedingF;

/// Score **one sample** against every candidate genotype, then **replace** that sample's
/// expected allele copies with the answer — one sample's share of the E-step.
///
/// **It is not idempotent, and the destructive input is the one thing a call site cannot
/// see.** The sample's old expected copies are read first, as the leave-one-out term of step
/// 1, and overwritten in step 4; calling it twice without an intervening M-step scores the
/// sample against a summary its own first call already moved.
///
/// Four things happen, in this order, and each is the next one's input:
///
/// 1. **This sample's concentration**, from the locus's seed plus what the *other* samples
///    showed here. The sample's own expected copies are subtracted from the cohort's total
///    first, so its own reads cannot reach it twice — once through its likelihood row and
///    once through the frequency they helped set
///    (`doc/devel/ng/spec/calling_priors.md` §6).
/// 2. **This sample's log-prior** over the candidate genotypes, from that concentration.
///    Which prior is a seam: two are shipped and they disagree by 11 points of genotype
///    accuracy on GIAB at 5×, each sample called on its own — 83.6% against 94.6%
///    (`doc/devel/ng/spec/calling_priors.md` §2.2) — so the model arrives as an argument
///    rather than being named here.
/// 3. **This sample's posterior** over the candidate genotypes: likelihood plus prior, then
///    normalised so the row sums to one.
/// 4. **This sample's expected allele copies**, which is the posterior read as copies —
///    every genotype's copies of each allele, weighted by how probable that genotype is.
///    Fractional, and never a call: that is what lets the loop work where no sample's
///    genotype is certain (`doc/devel/ng/spec/calling_em_loop.md` §1.3).
///
/// **No line of this function branches on the cohort size, and step 1 is why.** At one
/// sample the subtraction there is between a number and itself, so the concentration comes
/// back as the seed bit for bit and the loop reaches its fixed point on the first pass — by
/// arithmetic, not by a test of the cohort size
/// (`doc/devel/ng/spec/calling_em_loop.md` §7).
///
/// **The normalisation subtracts the largest score before exponentiating**, which is what
/// makes it safe at both ends: no term can overflow, because every exponent is at or below
/// zero; and the total is at least one, because the largest score's own term is exactly
/// `exp(0)`, so the division that follows cannot divide by zero. The row therefore always
/// has a maximum to subtract, which the genotype prior guarantees by flooring an impossible
/// genotype near −691 rather than writing `−∞`
/// (`doc/devel/ng/spec/calling_priors.md`; `calling_em_loop.md` §8).
///
/// **Nothing allocates.**
///
/// # Panics
///
/// Every check below is **held in release**, because this module's design makes a caller bug
/// an assertion rather than a `Result` (`doc/devel/ng/spec/calling_em_loop.md` §8). They do
/// not all earn that the same way, and the difference is worth having in the text, because a
/// later step trimming "redundant" checks has no other way to tell them apart:
///
/// - **the likelihood row and the posterior row must be one entry per candidate genotype,
///   and the copy table one row per genotype.** These three fail by *truncation* rather than
///   by crashing — `zip` and `chunks_exact` stop at the shorter side, so the genotypes past
///   the end keep whatever the previous sample left there and the call comes out confident
///   and wrong. Measured with the posterior-row check removed: a three-genotype locus scored
///   through a two-entry row returns copies `[1.333…, 0.667…]` where the answer is
///   `[1.0, 1.0]`, and **both** invariants the three-sample test asserts are still satisfied.
/// - **the seed must be one entry per allele.** This one is a message, not a catch: with it
///   removed, every mis-shaped seed is still refused in release a few lines later, by
///   [`fill_sample_concentration`] or by [`PriorRow::new`]. It is here so the panic names the
///   seed rather than the buffer that disagrees with it.
/// - **the sample's own expected copies must all be finite and non-negative.** Its own
///   sentinel makes this necessary: [`prepare_for_locus`] fills that row with
///   `UNWRITTEN_SCRATCH_VALUE`, which is `NaN`, and the leave-one-out term's `max(0, ·)`
///   *returns the other operand* on a `NaN` — so in release a row nobody wrote is absorbed
///   as a zero cohort term and the sample is scored against the bare seed, with the cohort's
///   evidence silently absent. The check that would have caught it downstream
///   (`SampleAlleleCopies::new`'s) is debug-only.
/// - **the total weight of the normalised row must be at least one.** The largest score's own
///   term is `exp(0)`, so it cannot be less — unless a score is `NaN`, which is what this
///   catches, at one comparison per sample per pass.
/// - **the largest score must be finite.** Note what this does *not* do: `largest_score` is
///   only ever assigned through `score > largest_score`, and every comparison against a `NaN`
///   is false, so **the maximum is never itself a `NaN`**. This check sees `±∞` — a row that
///   is entirely `NaN` reaches it as the `−∞` it started from. A `NaN` in any *other*
///   genotype is caught by the total-weight check above and by nothing else.
///
/// The remaining shapes are checked by [`fill_sample_concentration`] and [`PriorRow::new`],
/// which is why they are not repeated here.
// **`pub(crate)` with the dead-code lint expected, rather than `pub`.** Step D1 of
// `doc/devel/ng/impl_plan/calling_loop.md` is the caller; until it lands the only callers are
// this module's tests, and widening the visibility to silence the lint would be an accident
// dressed as a decision. `expect` rather than `allow`, and only outside the test build where
// the lint genuinely fires, so that D1 adding a real caller turns this line into a compile
// error and whoever writes it deletes the expectation.
#[cfg_attr(
    not(test),
    expect(dead_code, reason = "step D1 of the calling-loop plan is the caller")
)]
pub(crate) fn score_one_sample(
    buffers: SampleScoringBuffers<'_>,
    genotypes: &GenotypeTableView<'_>,
    prior: &dyn GenotypePriorModel,
    inbreeding: InbreedingF,
) {
    let SampleScoringBuffers {
        sample,
        seed_concentration,
        cohort_expected_copies,
        genotype_likelihoods,
        sample_concentration,
        prior_per_allele_workspace,
        prior_row,
        posterior_row,
        sample_expected_copies,
    } = buffers;

    let genotype_count = genotypes.genotype_count();
    let allele_count = genotypes.allele_count();
    assert_eq!(
        genotype_likelihoods.len(),
        genotype_count,
        "one genotype likelihood per candidate genotype: the table holds {genotype_count} \
         genotypes and sample {sample}'s likelihood row holds {}",
        genotype_likelihoods.len()
    );
    assert_eq!(
        posterior_row.len(),
        genotype_count,
        "one posterior entry per candidate genotype: the table holds {genotype_count} \
         genotypes and the posterior row holds {}, scoring sample {sample}",
        posterior_row.len()
    );
    assert_eq!(
        seed_concentration.len(),
        allele_count,
        "one seed concentration per allele: the table is built over {allele_count} alleles \
         and the locus's seed holds {}",
        seed_concentration.len()
    );
    // **Debug-only, unlike the release-held checks around it, and the reason is that it
    // cannot fire today.** Step 4 reads this table with `chunks_exact`, so a short one would truncate the
    // fold silently — but every `GenotypeTableView` comes from `GenotypeTable::build`, which
    // asserts this same identity as it builds, so no caller can present an inconsistent one.
    // It is here for a view built some other way, and held in debug because holding it in
    // release would be a check the suite cannot reach and therefore cannot keep honest.
    debug_assert_eq!(
        genotypes.genotype_allele_counts().len(),
        genotype_count * allele_count,
        "the copy table is one row of {allele_count} alleles per genotype, so \
         {genotype_count} genotypes need {} entries and it holds {}",
        genotype_count * allele_count,
        genotypes.genotype_allele_counts().len()
    );
    // The one check on a buffer's *contents* rather than its length, and the sentinel is why:
    // `prepare_for_locus` leaves this row `NaN`, and the leave-one-out `max(0, ·)` below
    // returns the other operand on a `NaN` — so in release an unwritten row becomes a zero
    // cohort term and the sample is quietly scored against the bare seed.
    assert!(
        sample_expected_copies
            .iter()
            .all(|copies| copies.is_finite() && *copies >= 0.0),
        "sample {sample}'s own expected allele copies are counts of genome copies, so every \
         entry must be finite and at or above zero; the likeliest cause is that a pass \
         reached this sample before anything wrote them: {sample_expected_copies:?}"
    );

    // 1. The seed, plus what the other samples showed. The sample's own copies are read
    //    here and overwritten at the end, so the read has to come first.
    fill_sample_concentration(
        Concentration::new(seed_concentration),
        CohortAlleleCopies::new(cohort_expected_copies),
        SampleAlleleCopies::new(sample_expected_copies),
        sample_concentration,
    );

    // 2. The prior over genotypes that concentration implies. `PriorRow::new` checks the
    //    genotype table's three flat views against it, so a mis-shaped table is refused
    //    before any implementation is entered.
    {
        let mut row = PriorRow::new(
            Concentration::new(sample_concentration),
            genotypes.genotype_allele_counts(),
            genotypes.log_multinomial_coeffs(),
            genotypes.homozygous_alleles(),
            prior_per_allele_workspace,
            prior_row,
        );
        prior.fill_genotype_log_priors(&mut row, inbreeding);
    }

    // 3. Likelihood plus prior, normalised into a posterior over genotypes.
    let mut largest_score = f64::NEG_INFINITY;
    for ((slot, likelihood), prior_of_genotype) in posterior_row
        .iter_mut()
        .zip(genotype_likelihoods.iter())
        .zip(prior_row.iter())
    {
        let score = likelihood.get() + prior_of_genotype.get();
        *slot = score;
        if score > largest_score {
            largest_score = score;
        }
    }
    assert!(
        largest_score.is_finite(),
        "the largest of sample {sample}'s {genotype_count} genotype scores came out \
         {largest_score} — an infinity, or every score was NaN, since a NaN never wins the \
         comparison that picks the largest; a likelihood row or a prior row reached here \
         already wrong"
    );

    let mut total_weight = 0.0;
    for slot in posterior_row.iter_mut() {
        *slot = (*slot - largest_score).exp();
        total_weight += *slot;
    }
    // **Release-held, and it is the module's only `NaN` detector.** The finiteness check
    // above cannot be one: `largest_score` is assigned only through `score > largest_score`,
    // which every `NaN` loses, so a `NaN` in a genotype that is not the most probable one
    // never reaches it. This check does, because `exp(NaN)` is `NaN`, `NaN + x` is `NaN`, and
    // `NaN >= 1.0` is false. It costs one comparison per sample per pass, not per genotype.
    assert!(
        total_weight >= 1.0 && total_weight.is_finite(),
        "the largest score's own term is exp(0) = 1, so sample {sample}'s total weight \
         cannot come out below one: got {total_weight} over {genotype_count} genotypes, \
         which means a score was NaN"
    );
    for slot in posterior_row.iter_mut() {
        *slot /= total_weight;
    }

    // 4. The posterior read as expected allele copies. Genotypes in the table's own order,
    //    so the same evidence gives the same sum at any thread count.
    sample_expected_copies.fill(0.0);
    for (copies_per_allele, &genotype_probability) in genotypes
        .genotype_allele_counts()
        .chunks_exact(allele_count)
        .zip(posterior_row.iter())
    {
        for (slot, &copies) in sample_expected_copies.iter_mut().zip(copies_per_allele) {
            // Skipping the zeros is worth having and the reason is that a diploid genotype
            // carries copies of at most two alleles, so most of a wide locus's copy table is
            // zero — 90 of the 126 entries at six alleles. Measured against a stub prior,
            // where this function's own arithmetic is what is being timed, removing the skip
            // cost between 8% and 29% across twelve runs at two and at six alleles. Against
            // the prior that actually ships it disappears: 0–4%, inside the run-to-run
            // drift, because the marginalized Dirichlet prior is about 250–280 ns of a
            // 380–405 ns call at six alleles — which is `spec/calling_em_loop.md` §2's point
            // that the prior, not the loop, carries the expensive function.
            if copies != 0 {
                *slot += genotype_probability * f64::from(copies);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::genotype_prior::MarginalizedDirichletPrior;
    use crate::ng::calling::{CallingScratch, CandidateAlleles, GenotypeTable};
    use crate::ng::locus_generation::LocusKind;
    use crate::ng::types::{LogProb, Ploidy};
    use std::sync::Arc;

    /// A genotype prior that writes numbers the test chose.
    ///
    /// **A stub is what makes the E-step's own arithmetic checkable by hand.** Both shipped
    /// priors go through `lgamma`, so a hand-computed expectation for one of them would be
    /// a transcription of a calculator rather than a check on this function; the real ones
    /// are used below on the properties they can be held to instead. The concentration this
    /// sees is not thrown away — it is the scratch's own buffer, so the test reads it back
    /// and checks it there.
    struct FixedLogPriors(Vec<f64>);

    impl GenotypePriorModel for FixedLogPriors {
        fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, _inbreeding: InbreedingF) {
            let (_workspace, out) = row.scratch_and_out();
            assert_eq!(
                out.len(),
                self.0.len(),
                "the fixture's row is the wrong width"
            );
            for (slot, &value) in out.iter_mut().zip(&self.0) {
                *slot = LogProb(value);
            }
        }

        fn name(&self) -> &'static str {
            "log-priors chosen by the test"
        }
    }

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a diploid")
    }

    fn outbred() -> InbreedingF {
        InbreedingF::try_new(0.0).expect("an outbred sample")
    }

    /// A SNP locus over `alternatives + 1` alleles, with its diploid genotype table.
    fn generic_locus(alternatives: usize) -> (CandidateAlleles, Arc<GenotypeTable>) {
        let mut alleles = CandidateAlleles::new(Box::from(b"A".as_slice()), LocusKind::Generic);
        for base in b"TCG".iter().take(alternatives) {
            alleles.admit(Box::from(&[*base][..]));
        }
        let table = GenotypeTable::build(diploid(), alleles.len());
        (alleles, table)
    }

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
            &FixedLogPriors(vec![0.0, 0.0, 2.0 * two]),
            outbred(),
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
        for (slot, value) in scratch.sample_genotype_likelihoods_mut(0).iter_mut().zip([
            -1000.0,
            -1000.0 + two,
            -1000.0,
        ]) {
            *slot = LogProb(value);
        }

        score_one_sample(
            scratch.sample_scoring_buffers_mut(0),
            &view,
            &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            outbred(),
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
            &MarginalizedDirichletPrior,
            outbred(),
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
                &MarginalizedDirichletPrior,
                outbred(),
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
            &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            outbred(),
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
            &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            outbred(),
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
            &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            outbred(),
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
            &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            outbred(),
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
            &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            outbred(),
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
            &FixedLogPriors(vec![0.0, 0.0, 0.0]),
            outbred(),
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
            &MarginalizedDirichletPrior,
            outbred(),
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
            &MarginalizedDirichletPrior,
            outbred(),
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
}
