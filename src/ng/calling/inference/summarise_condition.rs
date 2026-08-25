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
use crate::ng::calling::genotype_prior::{
    CohortAlleleCopies, Concentration, GenotypePriorModel, PriorRow, SampleAlleleCopies,
    fill_sample_concentration,
};
use crate::ng::calling::{CohortSummingBuffers, SampleScoringBuffers};
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

/// Add every sample's expected allele copies up into the cohort's — the **M-step** of one
/// pass.
///
/// The other half of the loop: the E-step turned the cohort's summary into each sample's
/// genotype probabilities, and this turns those back into the summary the next pass will use
/// (`doc/devel/ng/spec/calling_em_loop.md` §2). It is the whole of the M-step — there is no
/// second quantity, because the allele frequencies are the only thing that moves while the
/// frequency loop runs (§5's table).
///
/// # The order of the sum is this function's, and that is the point
///
/// **Floating-point addition is not associative, so "the sum over the samples" is not one
/// number until the order is fixed.** This walks the table in ascending sample index, so the
/// same evidence gives the same cohort copies at any worker count, on any machine, in any
/// build — which is the whole of `doc/devel/ng/spec/calling_em_loop.md` §8's determinism
/// contract. That is why the table arrives whole rather than a row at a time: a caller handed
/// one row per call decides the order, and a caller that parallelised over samples would
/// decide it differently on every run.
///
/// **How far apart the orders actually land, measured on this module's own fixture:** three
/// samples carrying `1.0`, `2⁻⁵³` and `2⁻⁵³` copies of an allele sum forward to exactly `1.0`,
/// because each tiny addend rounds away against the one — and backward to
/// `1.0000000000000002`, because the two tiny ones meet each other first and their total no
/// longer rounds away. One unit in the last place. **A difference that size does not move an
/// argmax** — not over the six genotypes this module's fixtures use, and not over the 21 of
/// the six-allele locus spec §2 works through. That is why the test for this compares the
/// summed copies *bitwise* rather than comparing genotypes: a test whose observable is the
/// call passes against a sum with no fixed order and proves nothing (spec §13 test 2).
///
/// # A sample the candidate step set aside
///
/// **Not handled here, and it must be before the loop runs on real evidence.** Spec §5.0 sets
/// aside a sample whose own reads earned an allele the cap then cut: it is scored against
/// nothing and *contributes nothing to this sum*, because its posterior would sit over the
/// wrong set of alleles and would pull the locus's frequencies toward the reference by exactly
/// the samples carrying the rarest alleles. `LocusEvidence::Generic` carries the flag
/// (`GenericLocusSample::genotype_must_be_missing`) and keeps such a sample's *index*, so this
/// table has a row for it.
///
/// **What happens today is that such a row is never written and this function refuses the
/// locus**, loudly, because `prepare_for_locus` left it holding the `NaN` sentinel and a `NaN`
/// propagates through every addition into the finiteness check below. That is the right
/// failure while the exclusion is unbuilt — a wrong cohort frequency would be silent — but it
/// is not the answer. Step D1 assembles the loop and owns the choice between skipping those
/// rows here and never giving them one.
///
/// # Panics
///
/// Four checks, all **held in release**, because a caller bug in this module is an assertion
/// rather than a `Result` (`doc/devel/ng/spec/calling_em_loop.md` §8):
///
/// - **the cohort row must name at least one allele**, and **the sample count must be at
///   least one**. These two look like belt-and-braces and are not: a sample count of zero
///   *satisfies* the size check below, because `0 == 0 × alleles`, and would hand back an
///   all-zero cohort row that reads like a locus where nobody carries anything.
/// - **the table must be `sample_count × alleles`**, with the product taken by `checked_mul`
///   so that a wrapped multiply cannot admit a shape nobody asked for.
/// - **every cohort entry must come out finite and non-negative.** The load-bearing one: it
///   is the only thing standing between an unwritten sample row and a cohort summary that
///   every sample's next prior is built from. See the comment on it for what it does not
///   catch.
///
/// Nothing allocates.
#[cfg_attr(
    not(test),
    expect(dead_code, reason = "step D1 of the calling-loop plan is the caller")
)]
pub(crate) fn sum_cohort_expected_copies(buffers: CohortSummingBuffers<'_>) {
    let CohortSummingBuffers {
        sample_count,
        per_sample_expected_copies,
        cohort_expected_copies,
    } = buffers;

    let allele_count = cohort_expected_copies.len();
    assert!(
        allele_count > 0,
        "every locus has a reference allele, so a cohort row of no alleles is a scratch that \
         was never prepared for this locus"
    );
    assert!(
        sample_count > 0,
        "a cohort has at least one sample, so a sum over none is a run whose sample order \
         went missing"
    );
    // **`checked_mul`, not `*`.** `prepare_for_locus` guards this very product the same way,
    // and for the same reason: a plain multiply wraps in release, and a wrapped product that
    // happens to equal the table's length would let the sum run over a shape nobody asked for
    // and come back looking like a summary. In debug it panics with "attempt to multiply with
    // overflow", which names the arithmetic rather than the locus.
    let expected_entries = sample_count.checked_mul(allele_count).unwrap_or_else(|| {
        panic!(
            "a cohort of {sample_count} samples over {allele_count} alleles needs a \
             per-sample copies table longer than a usize can index"
        )
    });
    assert_eq!(
        per_sample_expected_copies.len(),
        expected_entries,
        "the per-sample copies are samples × alleles, sample-major: {sample_count} samples \
         over {allele_count} alleles is {expected_entries} entries and the table holds {}",
        per_sample_expected_copies.len()
    );

    // Ascending sample index, and within each sample ascending allele. Every allele's
    // accumulator therefore receives its addends in ascending sample order, which is the
    // order the determinism contract names.
    //
    // **The first row is copied in rather than added to a zeroed buffer**, and the difference
    // is not the arithmetic — `0.0 + x` is exactly `x` for every value here, so the two are
    // bit-identical. It is that a `fill(0.0)` would overwrite the scratch's
    // `UNWRITTEN_SCRATCH_VALUE` sentinel *before* the walk, so a walk that summed no rows at
    // all would hand back a row of zeros that reads like a real summary. That cannot happen
    // today, because `sample_count` is checked above — but it is exactly the shape step D1
    // would introduce if it implemented spec §5.0 by skipping set-aside rows here, and a
    // locus where every sample was set aside would then report the cohort as carrying
    // nothing. Seeding from a row that must exist makes the empty walk unrepresentable.
    let mut rows = per_sample_expected_copies.chunks_exact(allele_count);
    let first = rows
        .next()
        .expect("the shape check above admits at least one sample's row");
    cohort_expected_copies.copy_from_slice(first);
    for own_copies in rows {
        for (total, &copies) in cohort_expected_copies.iter_mut().zip(own_copies) {
            *total += copies;
        }
    }

    // **Release-held, and it is what makes the `NaN` sentinel reach anybody.** A sample row
    // nobody wrote is `UNWRITTEN_SCRATCH_VALUE`, and `NaN` propagates through addition, so one
    // check over the alleles catches an unwritten row anywhere in the table — at a cost of
    // `alleles`, not `samples × alleles`. Without it the cohort's summary goes on to build
    // every sample's prior on the next pass.
    //
    // **On the first pass at a locus only, and the difference matters to whoever writes D1.**
    // The sentinel is written by `prepare_for_locus`, once per locus, and the per-sample rows
    // are deliberately *not* re-armed between passes, because `score_one_sample` reads each
    // sample's previous copies as its leave-one-out term. So from pass 2 a row this pass did
    // not write holds the **previous pass's finite, plausible value** and is summed as though
    // it were current — measured `[1.0, 3.0]` where pass 1 gave `[2.0, 2.0]` and one sample's
    // E-step did not run. Finite, non-negative, wrong, silent. A loop that can skip a sample
    // mid-run — which is what spec §5.0's exclusion would be — needs a per-pass written mask,
    // and that is D1's to build.
    //
    // **Complete for a non-finite input and not for a negative one**, and the asymmetry is
    // worth knowing rather than glossing: `NaN` and `±∞` survive addition, so a single bad
    // entry anywhere reaches the total, but two finite rows of `-1.0` and `3.0` sum to a
    // perfectly acceptable `2.0`. Checking every input instead would cost
    // `samples × alleles` on the loop's hot path to catch a state no producer can reach —
    // expected copies are a posterior-weighted sum of non-negative copy counts. The check is
    // on the output because that is where it is cheap, not because it is exhaustive.
    assert!(
        cohort_expected_copies
            .iter()
            .all(|total| total.is_finite() && *total >= 0.0),
        "the cohort's expected allele copies are counts of genome copies, so every entry must \
         be finite and at or above zero after summing {sample_count} samples; a NaN here is a \
         sample row the E-step never wrote, since NaN survives every addition: \
         {cohort_expected_copies:?}"
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

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
    /// Run whole passes over a cohort — the E-step for every sample in index order, then the
    /// M-step — with the samples presented in whatever order `order` names.
    ///
    /// Returns each sample's most probable genotype **after the last pass**, in the order it
    /// was presented in, and the cohort's summed expected copies.
    ///
    /// **`passes` is the whole reason this helper exists rather than a one-pass one.** On the
    /// first pass every sample is scored against a cohort row this function wrote before the
    /// loop, so no genotype is downstream of `sum_cohort_expected_copies` and a test reading
    /// winners off pass 1 passes against an M-step that computes nothing. Measured during
    /// B2's review: with the summing loop deleted, the one-pass version of the permutation
    /// test below stayed green while three other tests failed. From pass 2 the winners are
    /// scored against the row the previous pass's M-step produced, which is the regime the
    /// coupling under test actually occurs in.
    ///
    /// The posterior is one reused buffer, so each sample's winner is read off before the next
    /// is scored — which is also why step C3 of the plan has to take the genotype quality as it
    /// scores rather than afterwards.
    fn run_passes(
        passes: usize,
        likelihoods: &[Vec<f64>],
        starting_copies: &[Vec<f64>],
        order: &[usize],
        alleles: &CandidateAlleles,
        view: &GenotypeTableView<'_>,
        seed: &[f64],
    ) -> (Vec<usize>, Vec<f64>) {
        assert!(passes > 0, "a pass count of zero scores nothing");
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(order.len(), alleles, view);
        scratch.seed_concentration_mut().copy_from_slice(seed);

        let mut cohort = vec![0.0; view.allele_count()];
        for &sample in order {
            for (total, copies) in cohort.iter_mut().zip(&starting_copies[sample]) {
                *total += copies;
            }
        }
        scratch
            .cohort_expected_copies_mut()
            .copy_from_slice(&cohort);

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
        for _pass in 0..passes {
            winners.clear();
            for index in 0..order.len() {
                score_one_sample(
                    scratch.sample_scoring_buffers_mut(index),
                    view,
                    &MarginalizedDirichletPrior,
                    outbred(),
                );
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
            &likelihoods,
            &starting_copies,
            &[0, 1, 2],
            &alleles,
            &view,
            &seed,
        );
        let (called_rotated, cohort_rotated) = run_passes(
            2,
            &likelihoods,
            &starting_copies,
            &[2, 0, 1],
            &alleles,
            &view,
            &seed,
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
        assert_eq!(half_an_ulp_at_one, 2.0_f64.powi(-53));
        assert_eq!(one_ulp_at_one, 2.0_f64.powi(-52));

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
    /// **This is the identity `spec/calling_em_loop.md` §5.0 will break**: once a sample the
    /// candidate step set aside stops contributing, the total becomes the ploidy times the
    /// *callable* count. Whoever builds that in step D1 should expect to edit this test, and
    /// its failing is the point.
    #[test]
    fn the_cohort_carries_the_ploidy_for_every_sample_it_summed() {
        let (alleles, table) = generic_locus(2);
        let view = table.view();
        let (likelihoods, starting_copies, seed) = three_disagreeing_samples();

        let (_winners, cohort) = run_passes(
            1,
            &likelihoods,
            &starting_copies,
            &[0, 1, 2],
            &alleles,
            &view,
            &seed,
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
    /// sample's next prior is built from. **It is also how a sample set aside by
    /// `spec/calling_em_loop.md` §5.0 fails today** — loudly, which is the right failure while
    /// the exclusion is unbuilt, and not the answer, which step D1 owes.
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
    /// where nobody carries anything, which the finiteness check accepts. Measured during
    /// B2's review: with both `> 0` guards deleted the whole of `ng::calling` stays green.
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
    /// **This is not a bug in the M-step; it is the reason spec §5.0's exclusion needs a
    /// per-pass written mask**, which step D1 owes. If a later step makes skipping a sample
    /// mid-run impossible, this test should start failing and be deleted with a note saying
    /// why.
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
}
