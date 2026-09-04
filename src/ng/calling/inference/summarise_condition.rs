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
use crate::ng::calling::genotype_prior::seed_generic::{VariantClass, fill_locus_concentration};
use crate::ng::calling::genotype_prior::{
    CohortAlleleCopies, Concentration, GenotypePriorModel, PriorRow, SampleAlleleCopies,
    fill_sample_concentration, fill_ssr_seed,
};
use crate::ng::calling::likelihood::copy_shares;
use crate::ng::calling::likelihood::generic::{
    assemble_genotype_log_likelihood_row, fill_generic_emissions,
};
use crate::ng::calling::likelihood::ssr::genotype_log_likelihood_row;
use crate::ng::calling::likelihood::ssr_emission::SsrEmissionModel;
use crate::ng::calling::likelihood::{FrozenContamination, ReadGroupCalibrations};
use crate::ng::calling::quality::{
    ArtifactTestCounts, score_best_genotype, score_uncorrected_site_quality,
};
use crate::ng::calling::{
    CallingScratch, CandidateAlleles, CohortSummingBuffers, ContaminationMixture, ErrorSpreadTable,
    ExpectedAlleleCopies, FrozenParameters, GenericLocusSample, GenericSampleEvidence,
    GenotypeTable, LocusEvidence, LocusInference, ReadGroupParameters, RepeatTractProvenance,
    SampleGenotypeCall, SampleScoringBuffers, SsrSampleEvidence, fill_batch_allele_copies,
    fill_contaminant_allele_frequencies, fill_error_spreads,
};
use crate::ng::locus_generation::SsrDetail;
use crate::ng::parameter_estimation::Provenance;
use crate::ng::parameter_estimation::joint::ssr_fit::Slippage;
use crate::ng::parameter_estimation::ssr::RepeatCount;
use crate::ng::types::{AlleleId, Genotype, InbreedingF, LogProb, Ploidy};
use std::iter::repeat_n;
use std::num::NonZeroU32;

use super::repeat_tract_parameters::{TractPrior, TractScoringFits, tract_candidates};
use super::slippage_refit::{PooledSlipCounts, RefitEmissionCache, largest_movement, refit_cells};
use super::{LocusGenotyper, RunnableCallingLoopConfig};

/// Which prior one pass scores against — **a value, not a code path**.
///
/// **The first pass through a locus has no prior at all, and it cannot simply be handed a flat
/// one.** The leave-one-out prior is built from the cohort's expected allele copies, and those
/// are what a *previous* pass produces; on the first pass there is none, so the buffers holding
/// them are still the scratch's `NaN` sentinel. A caller that passed a flat
/// [`GenotypePriorModel`] would therefore still run step 1 and read them.
///
/// **What that costs depends on which buffer is unwritten, and both cases are measured.** With
/// the sample's own copies unwritten — the ordinary first pass — it **panics, in release as
/// well as debug**, on this function's own release-held check that those copies are finite. If
/// only the *cohort* row is unwritten, nothing panics in release and the concentration comes
/// back as **the seed exactly**, because the leave-one-out `max(0, ·)` returns the other
/// operand on a `NaN`; probed at a seed of `[1.0, 0.5]`, the concentration came back
/// `[1.0, 0.5]`. **The second case is the seeded first pass
/// `doc/devel/ng/spec/calling_em_loop.md` §3 exists to prevent**, arrived at by accident rather
/// than by choice. So the choice is a variant here rather than a model behind the seam.
///
/// **Why the first pass is flat and not seeded**, in one paragraph, because it is the whole of
/// §3: the seed says a locus is almost certainly invariant — about `1 − 3θ/2` of the prior mass
/// on the homozygous reference against roughly `θ` on a heterozygote, a pull of about 30 Phred
/// at `θ = 0.001`. Apply that on the first pass at a locus where the reads are thin and every
/// sample carrying the variant is scored homozygous reference; their expected copies of the
/// alternative come out near zero, so the cohort term is near zero on the second pass, so the
/// prior is still just the seed. **The loop converges, and it converges to no-variant, having
/// never let the reads speak.** GATK's allele-frequency calculation starts flat for the same
/// reason, and so does production's own E-step
/// (`src/var_calling/posterior_engine.rs`, `EmStepPhase::FirstIteration`).
///
/// **It runs at the start of every outer round, not only the very first.** For the first round
/// the reason is that the number does not exist yet; for the later ones it is different — a
/// round begins with new slippage numbers, and the expected copies it would inherit were
/// converged under the *old* ones, so carrying them in seeds the new round at the old round's
/// answer. That is the same self-reinforcing start, with the previous round standing in for the
/// seed (spec §3).
///
/// **Inherited from GATK and from production, and never measured here — soft.** The argument is
/// arithmetic about the prior's size, not a count of calls that moved. It should bite hardest
/// where the read likelihood is weakest against a 20-to-30 Phred prior, which is the tomato
/// panel's corner at three reads a position, and hardly at all at 300. Spec §12's question 7 is
/// the measurement.
#[derive(Clone, Copy)]
pub(crate) enum PassPrior<'a> {
    /// **No prior at all** — every candidate genotype equally likely, and no concentration
    /// built, so nothing reads the cohort's expected copies. The first pass of every round.
    Flat,
    /// The seed plus what the **other** samples showed here — every pass but the first.
    LeaveOneOut {
        /// Which of the two shipped priors this run scores against.
        model: &'a dyn GenotypePriorModel,
        /// This sample's inbreeding coefficient, frozen by the parameter pre-pass.
        inbreeding: InbreedingF,
    },
}

impl PassPrior<'_> {
    /// Whether this pass builds a concentration and reads the cohort's expected copies.
    ///
    /// **One of the two places the arms are told apart** — this predicate and the `match` in
    /// [`score_one_sample`] — kept as a named method so that what it gates reads as a
    /// property of the pass rather than as an inline `matches!`.
    #[inline]
    fn reads_the_cohort(&self) -> bool {
        matches!(self, Self::LeaveOneOut { .. })
    }
}

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
///    (`doc/devel/ng/spec/calling_priors.md` §6). **Skipped entirely on a
///    [`PassPrior::Flat`] pass**, which is what makes the flat pass possible at all: the
///    cohort's copies do not exist yet on the first pass through a locus.
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
// **`pub(crate)` rather than `pub`**: the chain from here up to `SummariseConditionLoop` is
// this module's, and the arm is the only thing outside it a run needs. Until D1 wrote that arm
// the whole chain was dead in a non-test build and carried a `cfg_attr(not(test),
// expect(dead_code))` for it — the expectation rather than an `allow`, so that the first real
// caller would turn the line into a compile error. It did, and the line is gone.
pub(crate) fn score_one_sample(
    buffers: SampleScoringBuffers<'_>,
    genotypes: &GenotypeTableView<'_>,
    prior: PassPrior<'_>,
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
    // **These two are checked here because the flat arm reaches neither of the functions that
    // used to check them.** `prior_row`'s width was `PriorRow::new`'s to check and
    // `sample_expected_copies`' was `fill_sample_concentration`'s, and a flat pass calls
    // neither — so a mis-shaped buffer made the seeded arm panic in release and made the flat
    // arm return a **wrong posterior in silence**. Measured on a 3-genotype locus given a
    // 2-entry prior row: `[0.199, 0.399, 0.402]` against the right `[0.25, 0.5, 0.25]`, the
    // tail entry being the stale value the buffer arrived holding, carried through the
    // normalisation because the score loop's `zip` stops at the shortest row while the
    // normalising loop divides all of them.
    assert_eq!(
        prior_row.len(),
        genotype_count,
        "one prior entry per candidate genotype: the table holds {genotype_count} genotypes \
         and the prior row holds {}, scoring sample {sample}",
        prior_row.len()
    );
    assert_eq!(
        sample_expected_copies.len(),
        allele_count,
        "one expected-copies entry per allele: the table is built over {allele_count} \
         alleles and sample {sample}'s copies row holds {}",
        sample_expected_copies.len()
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
    //
    // **Only on a seeded pass.** On a flat pass this row is *expected* to hold the sentinel —
    // that is the whole situation the flat pass exists for — and step 4 overwrites it without
    // ever reading it.
    assert!(
        !prior.reads_the_cohort()
            || sample_expected_copies
                .iter()
                .all(|copies| copies.is_finite() && *copies >= 0.0),
        "sample {sample}'s own expected allele copies are counts of genome copies, so every \
         entry must be finite and at or above zero; the likeliest cause is that a pass \
         reached this sample before anything wrote them: {sample_expected_copies:?}"
    );

    // 1 and 2. The prior over genotypes, and — on a seeded pass only — the concentration it
    //          is built from. The flat arm writes a zero row rather than branching inside
    //          step 3, so the hot loop below keeps one spelling; the cost is `genotypes`
    //          stores against `genotypes` calls to `exp`.
    match prior {
        PassPrior::Flat => prior_row.fill(LogProb(0.0)),
        PassPrior::LeaveOneOut { model, inbreeding } => {
            // The seed, plus what the other samples showed. The sample's own copies are read
            // here and overwritten at the end, so the read has to come first.
            fill_sample_concentration(
                Concentration::new(seed_concentration),
                CohortAlleleCopies::new(cohort_expected_copies),
                SampleAlleleCopies::new(sample_expected_copies),
                sample_concentration,
            );
            // `PriorRow::new` checks the genotype table's three flat views against the
            // concentration, so a mis-shaped table is refused before any implementation is
            // entered. **The flat arm does not go through it**, which is why the copy table's
            // own width is checked above rather than being left to this call.
            let mut row = PriorRow::new(
                Concentration::new(sample_concentration),
                genotypes.genotype_allele_counts(),
                genotypes.log_multinomial_coeffs(),
                genotypes.homozygous_alleles(),
                prior_per_allele_workspace,
                prior_row,
            );
            model.fill_genotype_log_priors(&mut row, inbreeding);
        }
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

/// Whether the cohort's expected allele copies have stopped moving — **the loop's whole
/// stopping rule**, and the one place the division that makes it portable is written.
///
/// Two rows of expected copies, the one the pass just ending was scored against and the one
/// its M-step produced. The loop stops when the largest change between them, **divided by
/// the number of chromosomes in the cohort**, falls below `threshold`
/// (`doc/devel/ng/spec/calling_em_loop.md` §6).
///
/// # Why expected copies and not the frequencies the locus reports
///
/// **Because this is the quantity that feeds back.** The M-step produces it and the next
/// E-step's leave-one-out prior is built from it, so its movement is what says the loop has
/// not settled. Production tests exactly this, and its comment records what testing the
/// *reported* estimate cost: that number is scaled by a pseudocount and does not feed back,
/// so a larger pseudocount damped the delta and stopped the loop earlier on an otherwise
/// identical trajectory ([`posterior_engine.rs:2702`](../../../../src/var_calling/posterior_engine.rs)).
///
/// # Why the division, which is the easy thing to drop
///
/// **Expected copies are a count and the threshold is a fraction.** At one diploid sample the
/// cohort carries 2 chromosomes and at a thousand it carries 2,000, so a movement of `1e-3`
/// raw copies is a frequency shift of 5 in 10,000 in the first case and 5 in 10,000,000 in
/// the second. A criterion written on raw counts therefore **tightens by the cohort size**,
/// silently, across exactly the range this caller commits to
/// (`doc/devel/ng/spec/design_principles.md` §0). Dividing puts the movement on the same
/// `[0, 1]` frequency scale as the numbers the locus reports, which is what lets one
/// inherited constant mean one thing from a single sample to several thousand.
///
/// # Why `all(… < threshold)` and not a maximum
///
/// **The two spellings are the same arithmetic on finite rows and they part on a row nobody
/// wrote.** `prepare_for_locus` fills *both* cohort rows with
/// [`UNWRITTEN_SCRATCH_VALUE`](crate::ng::calling::UNWRITTEN_SCRATCH_VALUE), which is `NaN`,
/// and [`advance_cohort_expected_copies`](crate::ng::calling::CallingScratch::advance_cohort_expected_copies)
/// leaves the previous row holding it until a pass has actually advanced. Every comparison
/// against a `NaN` is false, so:
///
/// - `fold(0.0, f64::max)`, `fold(f64::NEG_INFINITY, f64::max)` and a hand-written
///   `if d > largest` all **discard** the `NaN` — `f64::max` is documented to return the other
///   argument when one side is `NaN`, and a `>` comparison against one is false. The first and
///   third hand back their seed of `0.0` and the second hands back `−∞`; all three are below
///   any threshold, so the locus reports itself settled after one pass having compared against
///   nothing;
/// - `all(|d| d < threshold)` returns **false**, which is what
///   [`previous_cohort_expected_copies`](crate::ng::calling::CallingScratch::previous_cohort_expected_copies)'s
///   own documentation promises: *"every comparison against it is therefore false"*.
///
/// Three of the four natural spellings are wrong in the same direction, and the failure is
/// silent — a genotype, flagged converged, from a loop that ran one pass. That is measured and
/// pinned by `the_fold_spellings_of_the_delta_settle_where_this_one_does_not`.
///
/// **[`run_frequency_loop`] itself never hands over such a row**, and the guarantee is for its
/// successors rather than for it: the prior-free initialisation's M-step writes finite copies
/// before the first swap, so by the time this is first called the previous row is a real
/// estimate. What the spelling protects is any later caller that compares before a pass has
/// advanced — step C3's final pass, and step D1's outer rounds, which restart the
/// initialisation and could reorder the swap.
///
/// # Why `.abs()`, which two alleles cannot show
///
/// **At a biallelic locus the two movements are equal and opposite**, because the expected
/// copies sum to the cohort's chromosome total on every pass — so a signed comparison gives the
/// same verdict as an absolute one and the `.abs()` looks like decoration. From three alleles
/// on it is not: the reference allele can fall by more than the threshold while every
/// alternative rises by less than it, and a signed comparison then calls a moving locus
/// settled. Pinned by `cohort_expected_copies_have_settled_refuses_a_fall_larger_than_every_rise`,
/// on a row measured at `[−0.0015, +0.0008, +0.0007]` against the `1e-3` threshold.
///
/// # Panics
///
/// Four checks, all **held in release**, because a caller bug in this module is an assertion
/// rather than a `Result` (`doc/devel/ng/spec/calling_em_loop.md` §8):
///
/// - **the previous row must name at least one allele**, and **the two rows must be the same
///   length**. Two checks rather than one, so a test can tell which fired: `all` over nothing
///   is `true` and `zip` stops at the shorter side, so a pair of empty rows — the shape an
///   unprepared scratch would hand over — reports every locus settled on its first pass, and a
///   short row settles on the alleles it does not have. This is the same failure mode the
///   sentinel exists to catch, reached by the one door the sentinel cannot cover.
/// - **the cohort must hold at least one sample.** Zero samples is zero chromosomes, which
///   turns every non-zero movement into an infinity and a zero movement into a `NaN`; neither
///   is below the threshold, so no locus would ever settle and every one would spend the whole
///   cap. **The chromosome count itself needs no check**, and that is the reason this takes a
///   [`Ploidy`] and a sample count rather than the product: a `Ploidy` is at least one by
///   construction, so the only way to a bad product is a cohort of no samples, and an infinite
///   or `NaN` chromosome count is no longer expressible.
/// - **the threshold must be finite and above zero.**
///   [`CallingLoopConfig::validate`](super::CallingLoopConfig::validate) already refuses one
///   that is not, and this is the check that says so where the value is used — the loop is
///   reachable from tests that build a threshold by hand. `> 0.0` alone refuses zero and `NaN`;
///   the `is_finite()` half is there for an infinity, which would report every locus settled on
///   its first pass.
#[must_use]
fn cohort_expected_copies_have_settled(
    previous_expected_copies: &[f64],
    current_expected_copies: &[f64],
    ploidy: Ploidy,
    sample_count: usize,
    threshold: f64,
) -> bool {
    assert!(
        !previous_expected_copies.is_empty(),
        "every locus has at least its reference allele, so a cohort row of no alleles is a \
         scratch that was never prepared for this locus"
    );
    assert_eq!(
        previous_expected_copies.len(),
        current_expected_copies.len(),
        "the two cohort rows are one entry per allele of the same locus: the previous pass's \
         row holds {} entries and this pass's {}",
        previous_expected_copies.len(),
        current_expected_copies.len()
    );
    assert!(
        sample_count > 0,
        "a cohort has at least one sample, so a convergence test over none is a run whose \
         sample order went missing"
    );
    assert!(
        threshold.is_finite() && threshold > 0.0,
        "the convergence threshold is a fraction on the frequency scale, so {threshold} is \
         not one; CallingLoopConfig::validate refuses it where a run sets it"
    );
    // Chromosomes, not samples: it is what puts the movement below on the frequency scale.
    // `u8` widens to `f64` exactly at every value it can hold; `usize` does so up to 2^53, so
    // the product is exact for any cohort a machine could hold.
    let cohort_chromosomes = f64::from(ploidy.get()) * sample_count as f64;
    previous_expected_copies
        .iter()
        .zip(current_expected_copies)
        .all(|(before, after)| (after - before).abs() / cohort_chromosomes < threshold)
}

/// What one run of the frequency loop reports about *how* it stopped — the two fields that
/// travel into [`LocusInference`](crate::ng::calling::LocusInference) unchanged.
///
/// **Both are emitted, and neither is an error.** A locus that runs out of passes is called
/// anyway: production retired its non-convergence error precisely so that one hard site could
/// not kill a cohort run (`doc/devel/ng/spec/calling_em_loop.md` §6). What the caller owes is
/// that the flag reaches the output, because a genotype from a loop that did not settle is a
/// weaker claim than one from a loop that did and nothing downstream can tell them apart.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[must_use]
pub(crate) struct FrequencyLoopOutcome {
    /// How many seeded passes ran. **At least one, and the prior-free initialisation is not
    /// one of them** — it produces the estimate pass 1 is compared against, so a locus that
    /// settles immediately reports one pass rather than zero
    /// (`doc/devel/ng/spec/calling_em_loop.md` §3).
    pub(crate) passes: u32,
    /// Whether the copies stopped moving, or the loop simply reached `max_passes`.
    pub(crate) converged: bool,
}

/// Run the frequency loop at one locus: the prior-free pass, then seeded passes until the
/// cohort's expected copies stop moving or the cap is reached.
///
/// **The innermost of the design's three loops** — the only one that repeats at the shipped
/// configuration (`doc/devel/ng/spec/calling_em_loop.md` §2).
///
/// **What it does with the genotype-likelihood table depends on `reassembly`, and that is the
/// whole of the parameter.** With `None` — a run whose parameter fit found no contamination — it
/// reads the table the scratch already holds and never touches it, because the assembled row
/// reads no allele frequency and cannot move. With `Some`, it assembles the table again at the
/// head of every pass, because `q(o)` — the contaminating population's frequency for the allele
/// an observation shows — is the locus's own number and moves with the loop
/// (`doc/devel/ng/spec/read_likelihoods.md` §3.6). **The emissions underneath are not
/// recomputed either way**: they read no frequency, so what a pass costs is the assembly, which
/// is one multiply and one add inside a logarithm the row was taking anyway (§6.1).
///
/// **What it leaves behind is as much the point as what it returns.** On return the scratch
/// holds the converged cohort expected copies, each sample's own expected copies, and the
/// posterior row of *the last sample scored* — so a caller that wants every sample's genotype
/// must score them again. That final pass is step C3's, and it is a pass rather than a
/// read-back because the posterior row is one reused buffer.
///
/// # The order within a pass, and why the swap sits between the two halves
///
/// 1. the **E-step** scores every sample against the cohort row as it stands — which is what
///    the *previous* pass's M-step wrote;
/// 2. the swap makes that row the previous one and hands back a `NaN`-filled buffer for this
///    pass;
/// 3. the **M-step** fills it;
/// 4. the two rows are compared.
///
/// Swapping before the E-step instead would hand every sample a row of `NaN`s to be scored
/// against, and swapping after the comparison would compare a row against itself.
///
/// **The settled test comes before the cap test, and that order is the difference between
/// two claims.** A locus whose last allowed pass is the one that settles is *converged*, not
/// capped: §6 makes the flag a statement about the locus, and reporting the cap there would
/// understate every genotype at the site. Pinned by
/// `run_frequency_loop_reports_converged_when_the_last_allowed_pass_settles`.
///
/// # One sample, and a thousand
///
/// **No line of this function branches on the cohort size**, and none should be written
/// (`doc/devel/ng/spec/calling_em_loop.md` §7). At one sample the prior's leave-one-out
/// subtraction is between a number and itself, so the concentration comes back as the seed and
/// the loop reaches its fixed point by arithmetic: pass 2's copies equal pass 1's bit for bit
/// and it stops with `passes = 2`. The wasted second pass is spec §12's question 6 and is an
/// optimisation, not a correctness rule.
///
/// # Panics
///
/// Every check is **held in release** (`doc/devel/ng/spec/calling_em_loop.md` §8). This
/// function adds none of its own: the inbreeding coefficients now live on the scratch, one
/// per prepared row, and
/// [`CallingScratch::inbreeding_coefficient_by_row`](crate::ng::calling::CallingScratch)
/// refuses a scratch whose row map is not one entry per prepared row. It is the same check the
/// old slice argument carried, moved onto the scratch, where it fires at the **first read** of
/// the row map rather than in this function's preamble.
///
/// **The two directions fail differently, which is why the check covers both, and neither
/// fails the way this paragraph used to say** (measured with the check downgraded, C3b's
/// review). A map one entry **short** panics at the walk's last read — `index out of bounds:
/// the len is 2 but the index is 2` — which is loud, but names a slice where the reader needs
/// the cohort. A map one entry **long** is the silent one: the walk is over the prepared row
/// count, so the surplus entry is never read at all and a scratch claimed for a different
/// locus runs to completion. The rest of the checks belong to the two halves this drives and
/// to [`cohort_expected_copies_have_settled`].
///
/// **The ploidy is the genotype table's, and there is no second source of it.** The table was
/// built for a `(ploidy, allele count)` shape and `prepare_for_locus` already refuses a scratch
/// that disagrees with it about the alleles; taking a `Ploidy` argument beside it would add a
/// number nothing compares. Measured on the three-sample fixture with a diploid table: a
/// `ploidy` argument of 64 returned `passes: 2, converged: true` against the true
/// `passes: 4`, identically in debug and release, with nothing asserting — a too-large ploidy
/// loosens the threshold by the ratio and claims convergence it did not reach.
pub(crate) fn run_frequency_loop<SsrEmissionScratch>(
    scratch: &mut CallingScratch<SsrEmissionScratch>,
    genotypes: &GenotypeTableView<'_>,
    model: &dyn GenotypePriorModel,
    config: &RunnableCallingLoopConfig,
    reassembly: Option<&TableReassembly<'_>>,
) -> FrequencyLoopOutcome {
    let row_count = scratch.row_count();

    // The initialisation: one E-step on the reads alone, then the M-step, which together
    // mint the expected copies the first seeded pass's prior is built from. It is **not**
    // counted as a pass — `passes` counts the passes that had a prior — and it runs at the
    // start of every outer round, not only the very first (§3).
    for row in 0..row_count {
        score_one_sample(
            scratch.sample_scoring_buffers_mut(row),
            genotypes,
            PassPrior::Flat,
        );
    }
    sum_cohort_expected_copies(scratch.cohort_summing_buffers_mut());

    let max_passes = config.max_passes.get();
    let mut passes = 0;
    loop {
        // **Every row, in row order, and the row order is the fixed one spec §8 requires**:
        // the rows are the run's sample order with the uncallable samples' gaps closed up.
        // The coefficient is read one row at a time rather than walked as a slice, because it
        // lives on the same scratch as the buffers the scoring writes — and
        // `inbreeding_coefficient_by_row` refuses a row map that is not one entry per prepared row, which
        // is what makes this a walk over every row rather than over as many as happened to
        // arrive.
        // **The genotype likelihoods are assembled again here, and only where something is
        // contaminated.** `q(o)` is the frequency of an observation's allele among the samples
        // that ran beside this one, so it moves with the copies the last M-step produced
        // (`doc/devel/ng/spec/read_likelihoods.md` §3.6). The emissions underneath are not
        // recomputed — they read no frequency — so what this costs is the assembly.
        if let Some(reassembly) = reassembly {
            reassembly.assemble(
                genotypes,
                ContaminantFrequencies::TheLoopsOwnEstimate,
                scratch,
            );
        }
        for row in 0..row_count {
            let inbreeding = scratch.inbreeding_coefficient_by_row()[row];
            score_one_sample(
                scratch.sample_scoring_buffers_mut(row),
                genotypes,
                PassPrior::LeaveOneOut { model, inbreeding },
            );
        }
        scratch.advance_cohort_expected_copies();
        sum_cohort_expected_copies(scratch.cohort_summing_buffers_mut());
        passes += 1;

        if cohort_expected_copies_have_settled(
            scratch.previous_cohort_expected_copies(),
            scratch.cohort_expected_copies(),
            genotypes.ploidy(),
            row_count,
            config.convergence_threshold,
        ) {
            return FrequencyLoopOutcome {
                passes,
                converged: true,
            };
        }
        if passes >= max_passes {
            // **Emitted with the flag, never dropped and never an error.** One hard locus
            // must not kill a cohort run (§6).
            return FrequencyLoopOutcome {
                passes,
                converged: false,
            };
        }
    }
}

/// **What one pass needs to assemble the genotype-likelihood table again**, for a run where
/// something is contaminated.
///
/// It exists because [`run_frequency_loop`] otherwise knows nothing about the locus's evidence
/// or the run's parameters — it walks scratch rows — and with contamination on it has to, since
/// the table it reads is no longer a constant across the loop.
///
/// **`None` at the call site is the uncontaminated run**, and that is the whole of the branch:
/// with no fraction fitted the assembled row reads no frequency, so the table the driver
/// assembled once before the loop is the table every pass reads.
///
/// # What moving `q(o)` gives up, and it is not nothing
///
/// **A loop that re-estimates part of its own likelihood from its own posteriors is no longer
/// plain expectation-maximization**, and loses that algorithm's guarantee that the data
/// likelihood cannot fall between passes. Nothing here relied on the guarantee — spec §13's
/// fourth test already forbids asserting a monotone movement, and §6 stops the loop on the
/// cohort's expected copies settling rather than on a likelihood — so no check changes. What
/// changes is that the pass count is a quantity nobody has measured on real data with
/// contamination on. **What is measured is one synthetic locus**: three diploid samples over
/// three alleles, four reads of each allele at each, one library apiece, at a 5% fraction —
/// **7 passes against the same evidence's 4 with no fraction fitted**. §6's cap of 50 is
/// inherited and unmeasured, and this is the first thing in the loop that can push toward it.
pub(crate) struct TableReassembly<'a> {
    evidence: &'a LocusEvidence<'a>,
    parameters: &'a FrozenParameters<'a>,
}

impl<'a> TableReassembly<'a> {
    /// The evidence and parameters one locus's later assemblies read.
    pub(crate) fn of(
        evidence: &'a LocusEvidence<'a>,
        parameters: &'a FrozenParameters<'a>,
    ) -> Self {
        Self {
            evidence,
            parameters,
        }
    }

    fn assemble<SsrEmissionScratch>(
        &self,
        genotypes: &GenotypeTableView<'_>,
        frequencies: ContaminantFrequencies,
        scratch: &mut CallingScratch<SsrEmissionScratch>,
    ) {
        assemble_genotype_likelihood_table(
            self.evidence,
            self.parameters,
            genotypes,
            frequencies,
            scratch,
        );
    }
}

/// **Eight of the artifact summary's nine numbers while they are still being pooled** — the
/// ninth is the allele the other eight are pooled *for*, which is chosen before any of this
/// and is never summed.
///
/// **Seven of the eight are counts and are summed as integers**, which is not a
/// micro-optimisation but the removal of a question: `f64` addition is not associative, so a
/// cohort's read totals summed in another order could differ in the last bits and the two
/// binomial tests downstream would see different inputs at a different worker count. Whole
/// reads summed as `u64` cannot. The eighth is fractional by construction — a genotype
/// carrying one copy of two expects half a sample's reads — so it is an `f64` summed in the
/// run's fixed sample order, the same order and the same reason as the M-step's
/// (`doc/devel/ng/spec/calling_quality.md` §9).
#[derive(Default)]
struct PooledArtifactCounts {
    reference_reads: u64,
    reference_forward_reads: u64,
    reference_placed_left_reads: u64,
    alternative_reads: u64,
    alternative_forward_reads: u64,
    alternative_placed_left_reads: u64,
    total_reads: u64,
    genotype_expected_alternative_reads: f64,
}

impl PooledArtifactCounts {
    /// Add one **called** sample's reads to the pool, and its call's share of the expected
    /// alternative reads.
    ///
    /// **A sample the candidate step set aside never reaches here**, and the exclusion is
    /// one decision rather than two: it has no called genotype to derive an expectation
    /// from, so counting its reads in the observed total while it contributes nothing to
    /// the expectation would manufacture an apparent *excess* — and since only a deficit is
    /// penalised, that would quietly weaken the test rather than break it loudly
    /// (`doc/devel/ng/spec/calling_quality.md` §6.3).
    ///
    /// **The sample's depth is its reads on the locus's alleles, and nothing else.** Two
    /// kinds of read are outside it, both deliberately: a *partial* read showed no allele —
    /// it says the sample carries *at least* this, not what it carries — and a read whose
    /// allele candidate selection dropped reaches the view as pooled error mass with no
    /// count beside it (`GenericSampleEvidence::unmatched_q_sum`). Production's depth is the
    /// same quantity, a sum over the alleles its record carries
    /// ([`qual_refine.rs:92`](../../../../src/vcf/qual_refine.rs)).
    fn add_called_sample(
        &mut self,
        evidence: &GenericSampleEvidence<'_>,
        primary_alternative: AlleleId,
        primary_alternative_copies: u32,
        ploidy: Ploidy,
    ) {
        let mut depth = 0_u64;
        for observation in evidence.supported {
            let reads = u64::from(observation.num_reads);
            depth += reads;
            // **Two `if`s over one allele, and the reference arm cannot be the primary
            // alternative's**, because `LocusInference::new` refuses a summary that names
            // the reference as its alternative.
            if observation.allele.is_reference() {
                self.reference_reads += reads;
                self.reference_forward_reads += u64::from(observation.forward_reads);
                self.reference_placed_left_reads += u64::from(observation.placed_left_reads);
            } else if observation.allele == primary_alternative {
                self.alternative_reads += reads;
                self.alternative_forward_reads += u64::from(observation.forward_reads);
                self.alternative_placed_left_reads += u64::from(observation.placed_left_reads);
            }
        }
        self.total_reads += depth;
        // How many of this sample's reads the *call* expects to carry the alternative: a
        // heterozygote expects half its depth, a homozygote all of it, a homozygous
        // reference none. The expectation is read from the genotypes and never from the
        // fitted frequency, which adapts to an artifact and would excuse it (§6.2).
        self.genotype_expected_alternative_reads +=
            (f64::from(primary_alternative_copies) / f64::from(ploidy.get())) * depth as f64;
    }

    /// The pool as the two artifact tests read it.
    fn into_summary(self, primary_alternative: AlleleId) -> ArtifactTestCounts {
        ArtifactTestCounts {
            primary_alternative,
            reference_reads: self.reference_reads as f64,
            reference_forward_reads: self.reference_forward_reads as f64,
            reference_placed_left_reads: self.reference_placed_left_reads as f64,
            alternative_reads: self.alternative_reads as f64,
            alternative_forward_reads: self.alternative_forward_reads as f64,
            alternative_placed_left_reads: self.alternative_placed_left_reads as f64,
            total_reads: self.total_reads as f64,
            genotype_expected_alternative_reads: self.genotype_expected_alternative_reads,
        }
    }
}

/// **Which allele the artifact tests treat as *the* alternative**: the non-reference allele
/// the most reads reached, pooled over the samples the locus is called on — or `None` where
/// there is no such allele to name.
///
/// `pooled_reads` is the scratch buffer this fills and reads, one entry per allele of the
/// locus; it arrives holding the previous locus's totals and is zeroed here.
///
/// **`None` is two different situations and the caller need not tell them apart**: a locus
/// whose candidate table is the reference alone — 27.4% of built loci on the 63-accession
/// tomato panel and 27.3% on HG002 at 30×
/// ([`SelectionVerdict::Selected`](crate::ng::calling::allele_candidates::SelectionVerdict))
/// — and a locus with alternatives that no read reached. Both leave the two tests with
/// nothing to weigh, and production returns its baseline unchanged in exactly these two
/// cases ([`qual_refine.rs:79`](../../../../src/vcf/qual_refine.rs)).
///
/// **Ties go to the lowest allele id**, because the fold keeps the first *strict* maximum —
/// the same rule and the same reason as the genotype quality's argmax
/// ([`score_best_genotype`]): the allele table's order is fixed, so a run must not depend on
/// which of two equally supported alternatives a comparison happened to keep.
///
/// # Panics
///
/// **Held in release**, and both are about an observation rather than a buffer:
///
/// - **every observation must name an allele this locus holds.** The candidate id is
///   selection's mapping of the merge's own allele index, and a mapping applied against the
///   wrong table produces ids past the end. Without this the buffer index panics with a
///   message naming a slice, which sends the reader to the scratch rather than to the join.
/// - **an observation's forward-strand and placed-left counts cannot exceed its read
///   count.** They are shares of the same reads, and both artifact tests read them as
///   fractions: a fraction above one reaches the binomial tail as a probability above one
///   and comes back as a penalty rather than as a failure.
fn pool_reads_and_pick_primary_alternative(
    per_sample: &[GenericLocusSample<'_>],
    pooled_reads: &mut [u64],
) -> Option<AlleleId> {
    pooled_reads.fill(0);
    for (sample, locus_sample) in per_sample.iter().enumerate() {
        if !locus_sample.is_callable() {
            continue;
        }
        for observation in locus_sample.evidence.supported {
            let allele = usize::from(observation.allele.get());
            assert!(
                allele < pooled_reads.len(),
                "sample {sample}'s reads name allele {allele} and this locus is called over \
                 {} alleles, so the observation was mapped against a different allele table",
                pooled_reads.len()
            );
            assert!(
                observation.forward_reads <= observation.num_reads
                    && observation.placed_left_reads <= observation.num_reads,
                "sample {sample}'s {} reads of allele {allele} carry {} on the forward \
                 strand and {} placed left: both are shares of those same reads, and the \
                 artifact tests read them as fractions",
                observation.num_reads,
                observation.forward_reads,
                observation.placed_left_reads
            );
            pooled_reads[allele] += u64::from(observation.num_reads);
        }
    }
    // The first strict maximum over the alternatives, so ties go to the lowest id; allele 0
    // is the reference and is never a candidate for it.
    //
    // PANIC-FREE: `best_allele` indexes `pooled_reads`, whose length is the genotype table's
    // `allele_count`, and `GenotypeTable::build` refuses an `allele_count` above
    // `MAX_ALLELE_COUNT` — so the index is at most 65,535 and always fits a `u16`.
    let (best_allele, best_reads) = pooled_reads.iter().enumerate().skip(1).fold(
        (0_usize, 0_u64),
        |(best_allele, best_reads), (allele, &reads)| {
            if reads > best_reads {
                (allele, reads)
            } else {
                (best_allele, best_reads)
            }
        },
    );
    (best_reads > 0).then(|| AlleleId(u16::try_from(best_allele).expect("an allele id fits a u16")))
}

/// **Score every sample once more, take its best genotype, and mint what leaves the locus** —
/// the final pass (`doc/devel/ng/spec/calling_em_loop.md` §2's `finally`).
///
/// It runs after the frequency loop has stopped, against the frequencies it settled on, and
/// it is **a pass rather than a read-back** for one reason: `CallingScratch`'s posterior row
/// is a single genotype-length buffer that every sample in turn is scored into, so by the
/// time the loop returns only the last sample's posterior still exists
/// (`doc/devel/ng/spec/calling_quality.md` §3.1).
///
/// # Three things are computed here because here is the last moment they can be
///
/// - **each sample's genotype quality**, taken by [`score_best_genotype`] from the posterior
///   row *as that sample is scored*, in the same walk that picks its winner. Computing it
///   after the pass would need a `samples × genotypes` posterior table kept for the whole
///   locus — a second buffer the size of the largest one already allocated, to produce one
///   number per sample (§3.1).
/// - **the site quality before its artifact correction**, from the genotype likelihood table
///   the loop leaves behind. That table is per-worker scratch overwritten at the
///   next locus, and the fold over it is quadratic in cohort size, so computing it
///   downstream would both carry about half a megabyte per locus in flight at 3,000 samples
///   and put a quadratic computation on the run's one serial thread (§3.2).
/// - **the nine pooled counts the artifact correction consumes.** Eight are the evidence's,
///   which is released with the locus; the ninth needs the calls, which is why the whole
///   summary is built in this pass rather than at the input edge (§3.3).
///
/// **The correction itself is not here** — it is a few dozen operations on those nine
/// numbers plus the baseline, and it belongs to the first ordered output stage, where it can
/// be re-run at two settings over the same called stream without re-running the loop (§3.4).
///
/// # What the pass does not recompute
///
/// **The cohort's expected allele copies are the loop's, carried out unchanged.** Deriving
/// them downstream from the called genotypes gives a different number, because a call has
/// already thrown away the uncertainty these counts still carry, and site filtering and
/// emission both read them (`doc/devel/ng/spec/calling_em_loop.md` §9). This pass overwrites
/// each *sample's* own copies as it scores it — that is the E-step's fourth step — and never
/// re-sums the cohort's.
///
/// **So the genotypes are one E-step further on than the frequencies they are reported
/// beside**, which is exactly what production's own final pass does. At convergence the two
/// differ by less than the threshold that stopped the loop, by definition of having stopped.
///
/// # A sample the candidate step set aside
///
/// It is emitted as [`SampleGenotypeCall::Missing`], **scored against nothing and with no
/// quality beside it** — the enum is what makes that expressible, since a missing call and a
/// low-confidence one are different claims and emission must not conflate them
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.0). Such a sample is in neither the artifact
/// counts nor the expectation they are compared against (§6.3).
///
/// **And it has no scratch row at all**, which is what keeps every cohort-shaped quantity over
/// one cohort: the rows are the run's sample order with such a sample's gap closed up, so the
/// M-step sums the rows there are, the convergence delta divides by the chromosomes those rows
/// carry, and the site quality's count axis runs over the same samples — none of them having
/// to be told to skip anything. This pass walks the run's samples against a **row cursor**,
/// and the cursor is the whole of the map back (spec §5.0, §9; the choice was D1's, and B2 and
/// C3b both recorded it as open until then).
///
/// # One sample, and a thousand
///
/// **No line branches on cohort size.** At one sample the loop has already reached its fixed
/// point, so this pass reproduces its last E-step; at a thousand the same walk runs a
/// thousand times. What does grow is the site quality's fold, quadratically
/// (`doc/devel/ng/spec/calling_quality.md` §9).
///
/// # What it allocates, which is not nothing
///
/// **One `Vec` of calls, one owned `Genotype` per called sample, and one copy of the cohort's
/// expected copies** — the locus's *output*, which is owned by the caller and outlives the
/// scratch. The arithmetic itself allocates nothing: every buffer it reads and writes is the
/// worker's (`doc/devel/ng/spec/calling_em_loop.md` §8).
///
/// # Panics
///
/// Every check is **held in release** (`doc/devel/ng/spec/calling_em_loop.md` §8). This
/// function adds two, and both are the same class as the one
/// [`run_frequency_loop`] carries — a positional join between two per-sample lists:
///
/// - **the inbreeding coefficients must be one per sample of the run**, in both directions.
///   The walk is over the coefficients, so a short slice would call some samples and silently
///   leave the rest without a call at all.
/// - **the callable samples must be exactly the rows the scratch was prepared for.** The
///   cursor and the table are kept in step by construction rather than by a stored map, and
///   this is what says the construction held: a table filled for one set of samples and read
///   against another would score one sample's reads against another's likelihoods, which is a
///   wrong genotype rather than a crash.
///
/// The rest belong to the four functions this drives: [`score_one_sample`],
/// [`score_best_genotype`], [`score_uncorrected_site_quality`] and
/// [`LocusInference::new`], which is where the two new fields are checked against the allele
/// table.
#[allow(
    clippy::too_many_arguments,
    reason = "the seam this feeds (LocusGenotyper::call_locus) already takes four of these — \
              the evidence, the parameters, the candidates and the scratch — and the five \
              that remain are the prior model, the loop's outcome and the two warrants a \
              locus carries, all of which D1 has in hand when it calls this; a bundle \
              invented to satisfy the lint would be a type nothing else names"
)]
pub(crate) fn summarise_final_pass<SsrEmissionScratch>(
    scratch: &mut CallingScratch<SsrEmissionScratch>,
    genotypes: &GenotypeTableView<'_>,
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    model: &dyn GenotypePriorModel,
    candidates: CandidateAlleles,
    outcome: FrequencyLoopOutcome,
    weakest_provenance: Provenance,
    repeat_tract: Option<RepeatTractProvenance>,
) -> LocusInference {
    let inbreeding_by_sample = parameters.inbreeding_coefficient_by_sample();
    assert_eq!(
        inbreeding_by_sample.len(),
        evidence.sample_count(),
        "the inbreeding coefficients are one per sample of the run: the evidence at {} covers \
         {} samples and {} coefficients arrived",
        evidence.region(),
        evidence.sample_count(),
        inbreeding_by_sample.len()
    );

    // The SNP/indel per-sample rows, where there are any. They carry two things this pass
    // needs and the repeat-tract path has neither: which samples the candidate step set
    // aside, and the per-allele read counts the artifact summary pools. At a tract what goes
    // wrong is slippage rather than strand or read position, and that is already inside the
    // read likelihood — `doc/devel/ng/spec/calling_quality.md` §8 leaves a tract's quality to
    // a sibling document that is not written.
    let generic_samples = match evidence {
        LocusEvidence::Generic {
            region: _,
            per_sample,
        } => Some(*per_sample),
        LocusEvidence::Ssr { .. } => None,
    };
    // **The rows the artifact summary pools from, and the allele it weighs against the
    // reference — one value, because the second cannot exist without the first.** Held as
    // two `Option`s and re-paired at the use site, the combination that cannot arise —
    // an alternative allele with no rows behind it — would pool nothing and say nothing,
    // which is the silent-wrong-answer shape rather than the loud one.
    let artifact_pool = generic_samples.and_then(|per_sample| {
        pool_reads_and_pick_primary_alternative(per_sample, scratch.pooled_allele_reads_mut())
            .map(|primary_alternative| (per_sample, primary_alternative))
    });

    // **The genotype table's ploidy, and there is no second source of it.** The table was
    // built for a `(ploidy, allele count)` shape and every genotype's copies come out of it;
    // taking `FrozenParameters`' would add a number nothing compares (C2's ruling). It stays
    // a `Ploidy` all the way to the division it is the divisor of, so that no arithmetic
    // below can be handed a zero.
    let ploidy = genotypes.ploidy();
    let mut pooled = PooledArtifactCounts::default();
    let mut calls = Vec::with_capacity(evidence.sample_count());

    // **The run's sample order against a row cursor.** The scratch holds one row per sample
    // the locus is called on, in this same order with the uncallable samples' gaps closed up
    // (spec §5.0), so the cursor *is* the map from a row back to the sample that filled it —
    // there is nothing else to keep in step, and a sample with no row never reaches the
    // scoring at all.
    let mut row = 0_usize;
    for (run_sample, &inbreeding) in inbreeding_by_sample.iter().enumerate() {
        // **The map the table build read, not the predicate it was built from.** Deriving
        // "which samples have rows" a second way and comparing only the two counts would let a
        // *permutation* through — the table filled for one sample and read for another — which
        // is a wrong genotype rather than a crash. Asking the map whose row this is joins them
        // at every sample, and asserting that against the candidate step's own ruling catches
        // the two coming apart in either direction.
        let this_sample_owns_the_next_row = scratch
            .run_sample_of_each_row()
            .get(row)
            .is_some_and(|&claimed_by| claimed_by == run_sample);
        let callable = is_callable(evidence, run_sample);
        assert_eq!(
            this_sample_owns_the_next_row,
            callable,
            "sample {run_sample} at {} is {} by the candidate step and {} row {row} of the \
             likelihood table, so the table this pass reads was filled for a different set of \
             samples",
            evidence.region(),
            if callable { "callable" } else { "ruled out" },
            if this_sample_owns_the_next_row {
                "holds"
            } else {
                "does not hold"
            }
        );
        if !callable {
            calls.push(SampleGenotypeCall::Missing);
            continue;
        }
        score_one_sample(
            scratch.sample_scoring_buffers_mut(row),
            genotypes,
            PassPrior::LeaveOneOut { model, inbreeding },
        );
        let (winner, genotype_quality) = score_best_genotype(scratch.posterior_row());
        // PANIC-FREE: `score_best_genotype` returns an index into the posterior row, which
        // `prepare_for_locus` sized from this same table's genotype count.
        let copies_of_each_allele = genotypes
            .allele_counts_of(winner)
            .expect("the winner is an index into the row this table's own width sized");
        if let Some((locus_samples, primary)) = artifact_pool {
            pooled.add_called_sample(
                &locus_samples[run_sample].evidence,
                primary,
                copies_of_each_allele[usize::from(primary.get())],
                ploidy,
            );
        }
        calls.push(SampleGenotypeCall::Called {
            genotype: mint_genotype(copies_of_each_allele),
            genotype_quality,
            // **Asked of the likelihood row, here, while this pass still has it.** The row is
            // per-sample scratch the next locus overwrites, so nothing downstream can recover
            // what the reads said as against what the prior assumed — which is the whole of
            // why the flag is minted at the point of scoring (`vcf::assemble`'s module note,
            // and `doc/devel/ng/spec/vcf_output.md` §7.1).
            reads_were_uninformative: reads_said_nothing(scratch.sample_genotype_likelihoods(row)),
        });
        row += 1;
    }
    assert_eq!(
        row,
        scratch.row_count(),
        "the evidence at {} names {row} callable samples and the scratch was prepared for {} \
         rows, so the table this pass read was filled for a different set of them. The \
         per-sample check above catches every disagreement but this one: rows claimed past the \
         end of the run are never asked whose they are",
        evidence.region(),
        scratch.row_count()
    );

    // **After the samples, and it could as well have been before them**: the fold reads the
    // genotype likelihood table, which the loop built and this pass only reads
    // (`doc/devel/ng/spec/calling_quality.md` §3.2).
    //
    // **⚠ At a repeat tract this is the one place an ordinary-site quantity is read**, and it
    // is worth naming because everything else on the tract path was moved off one. The prior
    // seed here is the run's *frequency* spectrum — how allele frequencies are spread across
    // the population — and the tract's own *length* spectrum has no bearing on it. That is not
    // a defect introduced by the tract branch: `calling_quality.md` §8 leaves a tract's site
    // quality to a sibling document that is not written, so a tract shares the SNP/indel fold
    // until that lands. It is stated here rather than left for a reader to discover, because
    // `population_diversity.md` §1's fourth goal names exactly this shape of failure.
    let site_quality = score_uncorrected_site_quality(
        scratch.site_quality_buffers_mut(),
        genotypes,
        parameters.prior_seed(),
    );

    let cohort_expected_copies =
        ExpectedAlleleCopies::new(scratch.cohort_expected_copies().to_vec(), &candidates);
    LocusInference::new(
        evidence.region(),
        candidates,
        calls,
        cohort_expected_copies,
        outcome.converged,
        outcome.passes,
        weakest_provenance,
        repeat_tract,
        site_quality,
        artifact_pool.map(|(_, primary)| pooled.into_summary(primary)),
    )
}

/// The owned genotype a row of the table's copy counts describes — allele `a` repeated as
/// many times as that genotype carries it.
///
/// **The loop's currency is the [`GenotypeIdx`](crate::ng::calling::GenotypeIdx) and the
/// output's is this**, and the conversion happens once per called sample on the final pass
/// rather than anywhere in the loop: a `Genotype` owns a boxed slice, so minting one per
/// sample per *pass* would be an allocation on the hot path for a value every pass but the
/// last throws away (`doc/devel/ng/arch/calling_em_loop.md` §2).
fn mint_genotype(copies_of_each_allele: &[u32]) -> Genotype {
    let mut alleles = Vec::with_capacity(copies_of_each_allele.iter().sum::<u32>() as usize);
    for (allele, &copies) in copies_of_each_allele.iter().enumerate() {
        // PANIC-FREE: as at `pool_reads_and_pick_primary_alternative` — the row is one entry
        // per allele, and `GenotypeTable::build` caps `allele_count` at `MAX_ALLELE_COUNT`.
        let allele = AlleleId(u16::try_from(allele).expect("an allele id fits a u16"));
        alleles.extend(repeat_n(allele, copies as usize));
    }
    Genotype::new(alleles)
}

/// Which class of variant this locus is, for the prior's seed.
///
/// **Read off the candidate sequences and nowhere else.** A locus every one of whose
/// alternatives is the reference's own length is a substitution; anything else carries an
/// insertion or a deletion. The two classes take the same seed today — the projection reads
/// the shape of variation off allele counts the pre-pass does not separate by class — so this
/// moves no number yet. It is derived rather than defaulted so that the day the split arrives
/// it arrives at every locus at once (`doc/devel/ng/spec/calling_priors.md` §4.2, settled
/// 2026-08-22).
fn variant_class_of(candidates: &CandidateAlleles) -> VariantClass {
    let reference_length = candidates.reference().len();
    if candidates
        .iter()
        .all(|allele| allele.len() == reference_length)
    {
        VariantClass::Substitution
    } else {
        VariantClass::InsertionOrDeletion
    }
}

/// Whether the locus can be called on this run sample, or the candidate step ruled it out.
///
/// **A repeat tract rules no sample out** (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1): a
/// discovery round there can put back a length the cap cut, so no sample is locked out of the
/// locus for the rest of its calling.
fn is_callable(evidence: &LocusEvidence<'_>, run_sample: usize) -> bool {
    match evidence {
        LocusEvidence::Generic {
            region: _,
            per_sample,
        } => per_sample[run_sample].is_callable(),
        LocusEvidence::Ssr { .. } => true,
    }
}

/// **Whether a sample's reads said nothing about which genotype it has** — every genotype
/// equally probable under its own likelihoods, so the prior decides alone.
///
/// **Exact equality, not a tolerance, and the shape of the arithmetic is why.** A sample with
/// no reads at this locus contributes an empty sum to every genotype, and an empty sum is
/// `0.0` for all of them by the same code path — not three numbers that happen to agree. So a
/// tolerance would buy nothing here and would start calling *nearly* flat rows flat, which is
/// a different claim and one nothing has measured a threshold for.
///
/// **A row of one genotype is flat**, trivially and correctly: a haploid locus over one allele
/// has one genotype, the reads cannot distinguish it from anything, and there is nothing for
/// them to have said.
fn reads_said_nothing(likelihoods: &[LogProb]) -> bool {
    likelihoods.iter().all(|score| *score == likelihoods[0])
}

/// **The weakest warrant behind any parameter that reached this locus** — what the record is
/// entitled to claim about how well founded its genotypes are.
///
/// **Combined, never branched on** (`doc/devel/ng/spec/read_likelihoods.md` §4.4): a call
/// resting on one fitted parameter and one defaulted one is a defaulted call, and saying
/// otherwise launders the weaker of the two. [`Provenance::weaker_of`] is the ladder, and it
/// is the ladder `parameter_estimation` already states rather than one invented here.
///
/// **On the SNP/indel path, only the read groups whose reads are actually here.** A run's other
/// libraries contributed nothing to this locus, and charging it for a library that sent no read
/// would make a locus's warrant a property of the run rather than of the evidence. Both a
/// sample's whole-span observations and its partial ones name a read group, and both are
/// scored, so both count. **A repeat tract does not hold to that**, and the ⚠ below says why
/// and what it costs.
///
/// **What is not in it yet, and it is not nothing.** The prior's fitted spectrum carries no
/// provenance at all ([`SpectrumSeed`](crate::ng::calling::genotype_prior::SpectrumSeed) is
/// three numbers), on either path — a repeat tract's prior shape carries the rung it came from,
/// but that rung travels on the locus's own record rather than through this fold, because it is
/// a statement about the *prior* where every rung here is a statement about the *reads*.
///
/// # Each path folds the parameters its own row reads, and the two lists do not overlap
///
/// **A SNP/indel locus's row reads the calibrations** — one fitted error scale per read group —
/// and nothing else this fold covers.
///
/// **A repeat tract's row never reads a calibration at all.** It is scored from a stutter model
/// and a substitution rate per `(read group, candidate)`, which the tract's own parameter
/// assembly gathers and whose warrants it already folded
/// ([`TractScoringFits::weakest_warrant`](super::repeat_tract_parameters::TractScoringFits::weakest_warrant));
/// the calibration scale does not enter a tract's likelihood at all, so folding it in would
/// charge a tract for a parameter it never read.
///
/// **⚠ A tract's fold runs over every read group of the run, not the ones whose reads reached
/// it, and that is a departure from `doc/devel/ng/spec/read_likelihoods.md` §4.4** — which asks
/// for the weakest warrant of any parameter that *entered* the locus. The parameter table this
/// folds is built on the run's read-group axis (that type's own documentation says so and why),
/// so a tract in a run of many libraries can be reported `Defaulted` on account of a library
/// that sent it no read.
///
/// **It is the conservative direction** — a call is never claimed better founded than it is —
/// which is why it is stated rather than fixed here. Narrowing it means narrowing what the
/// parameter table covers, since the fold has nothing finer to fold; the read groups a tract's
/// observations name are what it would be narrowed to.
///
/// **A locus no read reached comes back [`Provenance::FittedHere`]** on the SNP/indel path,
/// because nothing weaker entered it: every sample is decided by the prior alone, and the prior
/// has no warrant to report. That is the fold's identity rather than a claim, and it is the one
/// answer here that would change the day the seed carries one. **At a tract it is not reached
/// at all**: the parameter table covers `read groups × candidates` whatever the reads did, so
/// a tract with no read still folds over a full table.
///
/// # Panics
///
/// At a repeat tract, if `tract_fits` holds no gathered tract — which says the parameters were
/// never read for this locus. On the SNP/indel path `tract_fits` is not read.
fn weakest_warrant_at_the_locus(
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    tract_fits: &TractScoringFits,
) -> Provenance {
    let calibration = parameters.calibration_by_read_group();
    let mut weakest = Provenance::FittedHere;
    match evidence {
        LocusEvidence::Generic {
            region: _,
            per_sample,
        } => {
            for locus_sample in per_sample.iter().filter(|sample| sample.is_callable()) {
                let read_groups = locus_sample
                    .evidence
                    .supported
                    .iter()
                    .map(|observation| observation.read_group)
                    .chain(
                        locus_sample
                            .evidence
                            .partials
                            .iter()
                            .map(|partial| partial.read_group),
                    );
                for read_group in read_groups {
                    let at = read_group.get() as usize;
                    let warrant = calibration
                        .get(at)
                        .unwrap_or_else(|| {
                            panic!(
                                "read group {at} has no calibration; the run supplied {}",
                                calibration.len()
                            )
                        })
                        .provenance;
                    weakest = weakest.weaker_of(warrant);
                }
            }
        }
        // **The tract's warrants are already folded, one cell per `(read group, candidate)`,
        // by the assembly that read them** — over the parameters a tract's row actually
        // scores from, which are not the calibrations the arm above reads.
        LocusEvidence::Ssr { .. } => {
            weakest = weakest.weaker_of(tract_fits.weakest_warrant());
        }
    }
    weakest
}

/// The SNP/indel path's per-sample evidence.
///
/// **Its callers have already chosen the path**, so a repeat tract reaching here is a branch
/// that was not taken rather than a case nobody wrote: the two paths' evidence is two different
/// shapes, and each of this driver's steps picks between them before it reads either.
fn generic_evidence_of<'a>(evidence: &'a LocusEvidence<'a>) -> &'a [GenericLocusSample<'a>] {
    match evidence {
        LocusEvidence::Generic {
            region: _,
            per_sample,
        } => per_sample,
        LocusEvidence::Ssr { region, .. } => unreachable!(
            "the repeat tract at {region} reached the SNP/indel path's evidence: its caller \
             branches on the evidence before reading it, so one branch is missing rather than \
             one path unbuilt"
        ),
    }
}

/// **Everything a repeat tract carries beside its samples' reads** — handed out together
/// because the two per-locus facts below are read within a line or two of the observations.
///
/// *(One caller today: the tract's table fill. The per-locus prior seed reads the same three
/// things and destructures the evidence itself, because it runs before the seed's own branch
/// has anything else to do with them.)*
struct TractEvidence<'a> {
    /// One entry per sample of the run, in run order. **A tract rules no sample out**
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1), so every sample has one.
    per_sample: &'a [SsrSampleEvidence<'a>],
    /// The tract's repeat unit and its two flanks.
    detail: &'a SsrDetail,
    /// How many whole repeats each candidate carries, parallel to the locus's candidate table.
    candidate_repeat_counts: &'a [NonZeroU32],
}

impl<'a> TractEvidence<'a> {
    /// **What the genotype prior believes about this tract's lengths, looked up once.**
    ///
    /// Both the prior's own seed and the third term of the read-likelihood mixture read it, and
    /// they must read the same one: the run reports which rung of the tract ladder answered, so
    /// a second lookup keyed differently would make the run's own record false. This is the one
    /// lookup.
    ///
    /// **The reference tract's repeat count is entry 0**, because the reference allele is id 0
    /// of every candidate table and these counts run parallel to it. It is what the spectrum's
    /// offsets are measured from, and passing a candidate's instead would re-centre the shape on
    /// that candidate.
    fn prior(&self, parameters: &'a FrozenParameters<'a>) -> TractPrior<'a> {
        let reference_repeats = RepeatCount(self.candidate_repeat_counts[0].get());
        TractPrior {
            reference_repeats,
            length_spectrum: parameters
                .ssr_length_spectrum_at(self.detail.motif.ssr_period(), reference_repeats),
        }
    }
}

/// The repeat-tract path's evidence, as [`generic_evidence_of`] is the SNP/indel path's — and
/// unreachable on that path for the same reason.
fn tract_evidence_of<'a>(evidence: &'a LocusEvidence<'a>) -> TractEvidence<'a> {
    match evidence {
        LocusEvidence::Ssr {
            region: _,
            per_sample,
            detail,
            candidate_repeat_counts,
        } => TractEvidence {
            per_sample,
            detail,
            candidate_repeat_counts,
        },
        LocusEvidence::Generic { region, .. } => unreachable!(
            "the SNP/indel locus at {region} reached the repeat-tract path's evidence: its \
             caller branches on the evidence before reading it, so one branch is missing rather \
             than one path unbuilt"
        ),
    }
}

/// **Everything at this locus that reads no allele frequency, computed once** — for every
/// claimed scratch row, and never again however many passes the loop takes.
///
/// **The name says what the two paths have in common**, because they fill different things:
/// what is shared is that no pass recomputes any of it
/// (`doc/devel/ng/spec/read_likelihoods.md` §6.1). What makes that checkable rather than merely
/// true is [`EmissionCost`](crate::ng::calling::EmissionCost), which this charges as it goes: a
/// version that recomputed the work on every pass would give identical genotypes, only slower,
/// so nothing but a counter can tell the two apart.
///
/// # What each path fills
///
/// **On the SNP/indel path, the emissions and nothing else.**
/// [`assemble_genotype_likelihood_table`] then folds them into a per-genotype row, once per
/// locus or once per pass depending on contamination. Three things are filled and none of them
/// reads a frequency: the locus's error-spread table, which is a property of the candidate
/// sequences and of the genotype being scored; each sample's charged error per observation; and
/// its compatibility verdict per `(partial read, candidate)`.
///
/// **At a repeat tract, the whole genotype-likelihood table**, because the tract's row builder
/// computes its emissions and assembles them in one call and there is no seam between the two
/// ([`genotype_log_likelihood_row`]). That costs nothing extra and changes no invariant: a
/// tract's row reads no allele frequency, so the table it writes is the table every pass reads
/// and [`assemble_genotype_likelihood_table`] has nothing left to do.
///
/// **That holds under contamination too, which is where the two paths part.** An ordinary
/// site's third mixture term is the cohort's own frequency for the allele an observation shows,
/// which the loop rewrites at every pass; a tract's is the *fit's* length spectrum for this
/// tract's stratum, frozen before calling starts
/// (`doc/devel/ng/spec/read_likelihoods.md` §4.5.1 against §3.6).
///
/// # Panics
///
/// Held in release (`doc/devel/ng/spec/calling_em_loop.md` §8). At a repeat tract, on a
/// candidate table that is not one tract's, on repeat counts that are not one per candidate,
/// and on a contamination list that disagrees with the fit the parameters were gathered under —
/// all three by
/// [`repeat_tract_parameters`](super::repeat_tract_parameters).
// **Test-only since the slippage round was built**: the driver now calls the `_for_round`
// spelling below, and this frozen wrapper is what the module's parity and cost fixtures
// exercise the two fills through. The `expect` rather than an `allow`, so that if the last
// test caller goes this line turns into a compile error instead of shipping dead code.
#[cfg_attr(not(test), expect(dead_code))]
fn fill_what_no_pass_recomputes<Model: SsrEmissionModel>(
    emission: &Model,
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    candidates: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    scratch: &mut CallingScratch<Model::Scratch>,
) {
    fill_what_no_pass_recomputes_for_round(
        emission,
        evidence,
        parameters,
        candidates,
        genotypes,
        scratch,
        SlippageRoundBuild {
            adopted_slippage: None,
            emissions: None,
        },
    );
}

/// [`fill_what_no_pass_recomputes`], with the slippage round's two channels — what the loop's
/// driver calls, since only it holds a round; everything else fills the frozen build above.
fn fill_what_no_pass_recomputes_for_round<Model: SsrEmissionModel>(
    emission: &Model,
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    candidates: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    scratch: &mut CallingScratch<Model::Scratch>,
    refit: SlippageRoundBuild<'_>,
) {
    match evidence {
        LocusEvidence::Generic { .. } => {
            // The slippage round is the repeat-tract path's: an ordinary site has no slippage
            // numbers to re-fit, so its rounds are ignored structurally (spec §5.1's closing
            // paragraph) and its build takes nothing from them.
            fill_generic_locus_emissions(evidence, parameters, candidates, genotypes, scratch);
        }
        LocusEvidence::Ssr { .. } => {
            fill_tract_likelihood_table(
                emission, evidence, parameters, candidates, genotypes, scratch, refit,
            );
        }
    }
}

/// **What one build of the tract's likelihood table takes from the slippage round, and hands
/// back to it** — nothing at the shipped configuration, where the round is frozen.
///
/// Two channels, each `None` on the frozen path:
///
/// - `adopted_slippage` — the round's re-fitted numbers, one entry per `(read group,
///   candidate)` cell of the gather, applied over the frozen gather so this round's table is
///   scored under them (spec §5.1's *rebuild the table*). `None` until a round has adopted
///   anything, which is also every build the frozen configuration ever makes.
/// - `emissions` — a cache the build copies each row's per-`(observation, candidate)`
///   emissions into as it fills them, so the round's attribution can read the same numbers
///   the genotype likelihoods were assembled from without recomputing an emission. `None`
///   wherever no round will run, and the build then copies nothing.
struct SlippageRoundBuild<'a> {
    adopted_slippage: Option<&'a [Option<Slippage>]>,
    emissions: Option<&'a mut RefitEmissionCache>,
}

/// The SNP/indel half of [`fill_what_no_pass_recomputes`]: the emissions, once per locus.
fn fill_generic_locus_emissions<SsrEmissionScratch>(
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    candidates: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    scratch: &mut CallingScratch<SsrEmissionScratch>,
) {
    let per_sample = generic_evidence_of(evidence);

    // **Once per locus rather than once per sample**: how far an allele's own error mass is
    // spread across the locus's others depends on the candidate sequences and on nothing a
    // sample showed.
    fill_error_spreads(candidates, genotypes, scratch.error_spreads_mut());

    scratch.charge_emission_build();
    let calibration = ReadGroupCalibrations::over(parameters.calibration_by_read_group());
    for row in 0..scratch.row_count() {
        let sample = per_sample[scratch.run_sample_of_each_row()[row]].evidence;
        scratch.charge_emission_row_fill(
            sample.supported.len() + sample.partials.len(),
            candidates.len(),
        );
        fill_generic_emissions(
            &sample,
            candidates,
            calibration,
            scratch.generic_row_mut(row),
        );
    }
}

/// **The repeat-tract half of [`fill_what_no_pass_recomputes`]: read this tract's fitted scoring
/// parameters, then score every sample's reads against every candidate genotype under them.**
///
/// One call per locus. It leaves the whole `rows × genotypes` genotype-likelihood table filled,
/// where the SNP/indel half leaves only the emissions underneath one.
///
/// # Why the parameters are read here rather than once per run
///
/// They are this tract's: which stutter model and which substitution rate a `(read group,
/// candidate)` cell takes depends on the candidate's own repeat count, so the table is gathered
/// per locus (`doc/devel/ng/spec/read_likelihoods.md` §4.4). It is gathered **inside** the
/// slippage round rather than before it, because that round exists to re-fit exactly these
/// numbers from the locus's own reads — a round that changed them and did not re-read them
/// would score the tract under the numbers it had just replaced.
///
/// # What it allocates, and it is the repeat path's own cost
///
/// **Two vectors per tract**: the candidate table as the row sees it — one
/// [`SsrCandidate`](crate::ng::calling::likelihood::ssr_emission::SsrCandidate) per candidate,
/// built straight from the evidence's repeat counts — and the scoring contexts,
/// `read groups × candidates` of them. **The contexts cannot live on the worker's scratch**:
/// they borrow the fits, which are a field of that same scratch, so a struct owning both would
/// refer to itself, and `#![forbid(unsafe_code)]` closes the usual escape.
///
/// A tract pays for a third vector elsewhere, in [`SummariseConditionLoop::call_locus`]: the
/// prior's seed builder takes the candidates' repeat counts as plain `u32`s where the evidence
/// carries them as non-zero ones, so the seed branch copies them. The SNP/indel path allocates
/// none of the three (`doc/devel/ng/spec/calling_em_loop.md` §8).
///
/// # Panics
///
/// On a candidate table that is not one repeat tract's, and on repeat counts that are not one
/// per candidate — both by
/// [`tract_candidates`](super::repeat_tract_parameters::tract_candidates).
///
/// On a contamination list that disagrees with the fit these parameters were gathered under —
/// by [`TractScoringFits::locus_parameters`](super::repeat_tract_parameters::TractScoringFits::locus_parameters).
fn fill_tract_likelihood_table<Model: SsrEmissionModel>(
    emission: &Model,
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    candidates: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    scratch: &mut CallingScratch<Model::Scratch>,
    mut refit: SlippageRoundBuild<'_>,
) {
    let tract = tract_evidence_of(evidence);
    let scored_candidates = tract_candidates(candidates, tract.candidate_repeat_counts);
    scratch.tract_fits_mut().gather_for_locus(
        &tract.detail.motif,
        &scored_candidates,
        tract.prior(parameters),
        parameters,
    );
    // **The slippage round's numbers go over the frozen gather, not instead of it.** The
    // gather is re-run each round because everything *else* a cell carries is frozen and this
    // is its one builder; what a round adopted then replaces the stutter models it re-fitted,
    // so this build scores the tract under the round's numbers (spec §5.1). `None` on every
    // frozen build.
    if let Some(adopted) = refit.adopted_slippage {
        scratch.tract_fits_mut().apply_refitted_slippage(adopted);
    }

    // **Charged before the buffers are taken, not as each row is scored.** The walk below runs
    // inside one borrow of this scratch, where no method of it can be called — so what a row
    // costs is charged here, from the same shapes the row is about to be handed.
    scratch.charge_emission_build();
    scratch.charge_table_assembly();
    let row_count = scratch.row_count();
    for row in 0..row_count {
        let run_sample = scratch.run_sample_of_each_row()[row];
        scratch.charge_emission_row_fill(
            tract.per_sample[run_sample].observations.len(),
            candidates.len(),
        );
        // **A tract's row assembles as it fills**, so the two halves are charged together —
        // one row assembled per row filled, where the SNP/indel path can charge them apart.
        scratch.charge_row_assembly();
    }

    // **The two walks below index the samples through the row map, and at a tract the map is
    // provably the identity** — `is_callable` is unconditionally true for a repeat tract, so
    // every sample claims a row in run order. Nothing here is therefore *tested* by a tract
    // whose rows and samples disagree, because no such tract exists; the map is read rather
    // than assumed so that the day one does, this walk needs no change.
    let contamination_of_each_read_group = parameters.contamination_by_read_group();
    let buffers = scratch.tract_locus_buffers_mut();
    let contexts = buffers.fits.scoring_contexts(&scored_candidates);
    let locus = buffers.fits.locus_parameters(
        &scored_candidates,
        &contexts,
        contamination_of_each_read_group,
    );
    let genotype_count = buffers.genotype_count;
    if let Some(cache) = refit.emissions.as_mut() {
        cache.begin_locus(candidates.len());
    }
    for row in 0..row_count {
        let run_sample = buffers.run_sample_of_each_row[row];
        let scored = row * genotype_count..(row + 1) * genotype_count;
        genotype_log_likelihood_row(
            emission,
            &tract.per_sample[run_sample],
            locus,
            genotypes,
            &mut buffers.genotype_likelihoods[scored],
            &mut *buffers.row_scratch,
        );
        // **Copied while this row's emissions are still in the one reused row scratch** — the
        // next row's fill overwrites them, and the slippage round's attribution reads the
        // cached emission per `(read, allele)` (spec §5.1). One `memcpy` per row, and only
        // where a round asked for it; the frozen path takes this branch never.
        if let Some(cache) = refit.emissions.as_mut() {
            cache.push_row(buffers.row_scratch.emissions());
        }
    }
}

/// **Attribute every read of one repeat tract to the candidate alleles, under the genotype
/// posteriors, and pool the slips** — the E-side of the slippage round's re-fit
/// (`doc/devel/ng/spec/calling_em_loop.md` §5.1).
///
/// A read's weight is split across the genotypes by each genotype's posterior, and within a
/// genotype across its alleles by the responsibility `copy share × cached emission`,
/// normalised over the genotype's own alleles. **Posteriors, not called genotypes** — the
/// spec's explicit ruling, which is HipSTR's choice and not production's hard assignment
/// (`src/ssr/cohort/em.rs:1192`): a call would push a thinly covered tract's numbers toward
/// whichever genotype won by a whisker.
///
/// The posteriors are minted here by scoring each row once against the converged frequencies
/// — the loop leaves only the last row's posterior behind, so a caller wanting every row's
/// must score them again (the same walk the final pass makes). **The scoring's one side
/// effect is undone before the next row**: `score_one_sample` replaces the row's expected
/// copies with the post-convergence posterior's, and the final pass will make the same walk
/// again later, so this one stashes each row's copies and puts them back — a round that
/// merely *measured* must leave the loop's state as it found it, or a re-fit whose numbers
/// never moved would still change the locus's calls.
///
/// # What is in the pool and what is out
///
/// - **A read that ran out inside the tract is out.** It witnessed a lower bound, not a
///   length, so there is no length difference to bin. Production's pool has no such reads
///   either — its per-sample length counts are over spanning reads.
/// - **A `(read, allele)` pair whose length difference is not a whole number of motif units
///   is out** — the part-repeat shares are placeholders the re-fit must not fit (spec §5.1,
///   *"hold the count at three"*). Out of numerator and denominator both: its weight joins
///   neither the slipped count nor the expected count.
/// - **A pair whose cell the fit has no frozen numbers for is out**, for the same shape of
///   reason: there is no frozen level to expect slips from and nothing to pull a re-fit
///   toward, and [`TractScoringFits::apply_refitted_slippage`] leaves such a cell's shipped
///   constants standing.
/// - **A read a genotype explains not at all** — every one of its alleles' emissions zero —
///   contributes nothing under that genotype: what explains it is the junk or contaminant
///   term, which attributes to no allele.
///
/// # Determinism
///
/// Rows in row order, observations in the merge's order within each, genotypes and alleles
/// in table order — every `f64` sum below runs in one fixed order, so the pooled counts are
/// identical at any worker count (spec §8).
fn attribute_reads_under_the_posteriors<SsrEmissionScratch>(
    scratch: &mut CallingScratch<SsrEmissionScratch>,
    genotypes: &GenotypeTableView<'_>,
    prior_model: &dyn GenotypePriorModel,
    tract: &TractEvidence<'_>,
    candidates: &CandidateAlleles,
    emissions: &RefitEmissionCache,
) -> PooledSlipCounts {
    let period = i64::from(tract.detail.motif.ssr_period().get());
    let candidate_count = candidates.len();
    let candidate_byte_lengths: Vec<i64> = candidates.iter().map(|b| b.len() as i64).collect();
    // Copied out because the walk below borrows the whole scratch mutably to score each row,
    // and the frozen cells live on that same scratch. One small vector per round, on the
    // round's own path only.
    let frozen_of_each_cell: Vec<Option<Slippage>> =
        scratch.tract_fits().frozen_slippage_cells().to_vec();
    let copy_share = copy_shares(genotypes.ploidy());
    let counts_of_each_genotype = genotypes.genotype_allele_counts();
    let allele_count = genotypes.allele_count();

    let mut pooled = PooledSlipCounts::default();
    let mut stashed_copies: Vec<f64> = Vec::new();
    for row in 0..scratch.row_count() {
        let inbreeding = scratch.inbreeding_coefficient_by_row()[row];
        stashed_copies.clear();
        stashed_copies.extend_from_slice(scratch.sample_expected_copies(row));
        score_one_sample(
            scratch.sample_scoring_buffers_mut(row),
            genotypes,
            PassPrior::LeaveOneOut {
                model: prior_model,
                inbreeding,
            },
        );
        let run_sample = scratch.run_sample_of_each_row()[row];
        let sample = &tract.per_sample[run_sample];
        let posterior = scratch.posterior_row();
        for (position, observation) in sample.complete_observations() {
            let reads = f64::from(observation.num_obs);
            let read_group = observation.read_group.0 as usize;
            for (genotype, &genotype_posterior) in posterior.iter().enumerate() {
                if genotype_posterior <= 0.0 {
                    continue;
                }
                let carried_copies =
                    &counts_of_each_genotype[genotype * allele_count..][..allele_count];
                // The genotype's own normaliser: what its alleles explain of this read,
                // copy-weighted — the same product the likelihood row sums
                // (`likelihood::ssr`'s row), read back from the build's cache.
                let mut explained = 0.0;
                for (candidate, &copies) in carried_copies.iter().enumerate() {
                    if copies == 0 {
                        continue;
                    }
                    explained += copy_share[copies as usize]
                        * emissions.emission_at(row, position, candidate);
                }
                if explained <= 0.0 {
                    continue;
                }
                for (candidate, &copies) in carried_copies.iter().enumerate() {
                    if copies == 0 {
                        continue;
                    }
                    let responsibility = copy_share[copies as usize]
                        * emissions.emission_at(row, position, candidate)
                        / explained;
                    let weight = reads * genotype_posterior * responsibility;
                    if weight <= 0.0 {
                        continue;
                    }
                    let cell = read_group * candidate_count + candidate;
                    let Some(frozen) = frozen_of_each_cell[cell] else {
                        continue;
                    };
                    let bp_diff =
                        observation.bases.len() as i64 - candidate_byte_lengths[candidate];
                    if bp_diff % period != 0 {
                        continue;
                    }
                    pooled.add(bp_diff / period, weight, frozen.level);
                }
            }
        }
        // Put the row's expected copies back — see the function note for why a measuring
        // round must not move the loop's state.
        scratch
            .sample_expected_copies_mut(row)
            .copy_from_slice(&stashed_copies);
    }
    pooled
}

/// Where `q(o)` — the contaminating population's frequency for the allele an observation shows —
/// comes from on one assembly of the table.
///
/// **Three answers and not two, because the loop has a state before its first pass.**
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
enum ContaminantFrequencies {
    /// The run's parameter fit identified no contamination anywhere, so there is no mixture and
    /// the row computes `doc/devel/ng/spec/read_likelihoods.md` §3.3 — the plain formula. This
    /// is what a single sample gets, and it is the *simple* case for that model rather than the
    /// weak one.
    NoneFitted,
    /// **The run fitted a fraction, but no pass has run**, so there is no frequency to read —
    /// and this is the table the prior-free initialisation pass scores against, whose one job is
    /// to turn the reads into a first estimate of the cohort's expected copies
    /// (`doc/devel/ng/spec/calling_em_loop.md` §3). It therefore scores the reads **alone**:
    /// §3.3's formula, which is what this model computes wherever `c` is zero.
    ///
    /// **The alternative that was tried first was a flat `q(o)`** — every candidate allele
    /// equally likely in the contaminating population — and it is worse for the reason §3
    /// rejects the seeded start: it does not say nothing, it damps the reads. `c · q` is a floor
    /// under every observation's mixture that no genotype can lower, so it compresses the
    /// differences between genotypes on the one pass whose whole purpose is to let those
    /// differences speak. Computed on §3.6's own formula at a hom-ref genotype scoring four
    /// alternative reads, `ε̄ = 0.01`, a spread of 3 and `c = 0.05`: a flat `q = 0.5` makes those
    /// reads **28 Phred cheaper** to explain than a converged `q = 0.05` does, against 3.7 Phred
    /// *dearer* for scoring them with no mixture at all. The reads-alone start is the closer of
    /// the two to where the loop settles, and it needs no rule the model does not already have.
    ///
    /// **What it decides is where the iteration starts, not what the model is** — every later
    /// assembly reads the loop's own estimate, contamination and all.
    TheReadsAlone,
    /// The loop's own current expected allele copies, summed over each sequencing batch and with
    /// the scored sample's own copies taken out of its own batch.
    TheLoopsOwnEstimate,
}

/// **Fold the whole `rows × genotypes` genotype-likelihood table out of the emissions**, one row
/// per claimed scratch row.
///
/// **Once per locus where nothing is contaminated, and once per pass where something is.** With
/// no fraction fitted the assembled row reads no frequency, so it is the same value at every
/// pass and the driver assembles it once; with a fraction fitted, `q(o)` moves with the loop and
/// the row has to be assembled again (`doc/devel/ng/spec/read_likelihoods.md` §3.6, §6.1).
///
/// # At one sample
///
/// **There is no mixture at all, and that is settled upstream rather than here.** Contamination
/// is a comparison between samples, so a single-sample run has nothing to fit a fraction from
/// and `RunParameters::view` routes it to `FrozenParameters::uncontaminated` — *emit it as
/// absent*, not a fitted zero, and the row then computes the plain formula. **This function
/// never sees that case as a contaminated one**, which is worth saying because the arithmetic
/// below would not refuse it: a lone sample is alone in its batch, the leave-one-out subtraction
/// empties that batch, and the no-evidence row would score it against the reference — close to
/// the plain formula but not equal to it, since an unexplained *reference* read would still pick
/// up `c` of the contaminant's mass. The named constructor is what keeps the two apart.
///
/// # At a repeat tract there is nothing here to do
///
/// **The tract's table is already filled**, by [`fill_what_no_pass_recomputes`], because the
/// tract's row builder computes its emissions and assembles them in one call — and that table
/// is final, because a tract's row reads no allele frequency. **Contaminated or not**: the
/// tract's third mixture term is the fit's length spectrum for its stratum, frozen before
/// calling starts, where an ordinary site's is the cohort's own estimate and moves with the
/// loop (spec §4.5.1 against §3.6). So a contaminated tract reaches this function and returns
/// from it having nothing to fold, and its table is assembled once like any other.
///
/// # Panics
///
/// Held in release (`doc/devel/ng/spec/calling_em_loop.md` §8), on a scratch whose contaminant
/// tables were not prepared for a run that fitted a fraction.
fn assemble_genotype_likelihood_table<SsrEmissionScratch>(
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    genotypes: &GenotypeTableView<'_>,
    frequencies: ContaminantFrequencies,
    scratch: &mut CallingScratch<SsrEmissionScratch>,
) {
    let allele_count = genotypes.allele_count();
    // **Held in release**, and the direction that makes it worth the check is the second one: a
    // contaminated run handed `NoneFitted` is scored on the plain formula at every locus of the
    // run, with the fraction the pre-pass fitted silently unused and nothing in the output
    // saying so. The other direction — an uncontaminated run asked for frequencies — cannot get
    // far, since its contaminant tables were never sized.
    assert_eq!(
        parameters.contamination_is_absent(),
        frequencies == ContaminantFrequencies::NoneFitted,
        "a run that fitted no contamination is the only one whose assemblies say so; a run that \
         did scores the reads alone on the initialisation assembly and its own estimate on every \
         later one. The caller picks the source from the parameters, so the two cannot disagree"
    );

    // **A repeat tract's table was filled whole by the emission build**, so there is nothing to
    // fold here — and that holds under contamination too, which is where the two paths part.
    //
    // **The SNP/indel mixture moves with the loop and a tract's does not.** `q(o)` at an
    // ordinary site is the cohort's frequency for the allele an observation shows, which the
    // loop rewrites at every pass (spec §3.6). A tract's third term is the *fit's* length
    // spectrum for this tract's stratum, frozen before calling starts — §4.5.1 chose it over
    // the cohort's own frequencies for exactly that reason, because contamination must not move
    // from one pass to the next. So an uncontaminated tract and a contaminated one are both
    // scored from one table, built once.
    let per_sample = match evidence {
        LocusEvidence::Generic {
            region: _,
            per_sample,
        } => *per_sample,
        LocusEvidence::Ssr { .. } => return,
    };

    match frequencies {
        // Neither reads a contaminant frequency, so neither touches the tables — and the
        // initialisation assembly leaves them holding the unwritten sentinel, which is what
        // stops it from being reached by a later assembly that forgot to fill them.
        ContaminantFrequencies::NoneFitted | ContaminantFrequencies::TheReadsAlone => {}
        // **Once per locus per pass**, and before any sample's own copies are taken back out of
        // it: what a batch holds is a property of the locus, and what one sample is scored
        // against is a property of that sample.
        ContaminantFrequencies::TheLoopsOwnEstimate => {
            let buffers = scratch.batch_copy_buffers_mut();
            fill_batch_allele_copies(
                buffers.expected_copies_by_run_sample,
                parameters.batch_of_each_sample(),
                buffers.allele_count,
                buffers.batch_allele_copies,
            );
        }
    }

    // **The run's half of the mixture is checked once here rather than once per row.** Those
    // checks are `read groups × batches`, and the row loop below runs them once a sample once a
    // pass — at a thousand libraries that is more work than the arithmetic the pass exists to
    // do. What is left per row is a check on the table that row's fill just wrote.
    let frozen_contamination =
        (frequencies == ContaminantFrequencies::TheLoopsOwnEstimate).then(|| {
            FrozenContamination::new(
                parameters.contamination_by_read_group(),
                parameters.batch_of_each_read_group(),
                parameters.batch_count(),
            )
        });
    let calibration = parameters.calibration_by_read_group();

    let row_count = scratch.row_count();
    scratch.charge_table_assembly();
    for row in 0..row_count {
        scratch.charge_row_assembly();
        let run_sample = scratch.run_sample_of_each_row()[row];
        let sample = per_sample[run_sample].evidence;
        if frequencies == ContaminantFrequencies::TheLoopsOwnEstimate {
            // **This sample leaves itself out of its own batch.** A contaminating read is
            // somebody else's by definition, so the population it is drawn against must not
            // include the individual being scored — the same subtraction the genotype prior
            // makes, one axis over.
            // **Asked by sample and not by index**, so that reaching for the read-group
            // batching here is a type error: `batch_of_read_group` takes a `ReadGroupId`. The
            // two batchings are the same slice type over different axes and the same length at
            // one library per sample, so a transposition otherwise passes every shape check.
            let own_batch = parameters.batch_of_sample(run_sample);
            let buffers = scratch.contaminant_frequency_buffers_mut(row);
            // **What comes back is how many batches had nothing left to read a frequency off**
            // — a locus scored against a contaminating population nobody measured. Reporting it
            // is E2b's, which is where the run says what contamination it used; it is bound
            // here rather than dropped so that the next reader sees it is owed rather than
            // handled.
            let _batches_with_no_evidence = fill_contaminant_allele_frequencies(
                buffers.batch_allele_copies,
                SampleAlleleCopies::new(buffers.own_expected_copies),
                own_batch,
                buffers.allele_count,
                buffers.contaminant_allele_frequencies,
            );
        }
        let buffers = scratch.generic_row_buffers_mut(row);
        let mixture = match frozen_contamination {
            None => ContaminationMixture::uncontaminated(),
            Some(frozen) => {
                frozen.with_frequencies(buffers.contaminant_allele_frequencies, allele_count)
            }
        };
        assemble_genotype_log_likelihood_row(
            &sample,
            genotypes,
            ReadGroupParameters::new(calibration, mixture),
            ErrorSpreadTable::over(buffers.error_spreads, genotypes),
            buffers.row_scratch,
            buffers.genotype_likelihoods,
        );
    }
}

/// **Summarise the cohort, then condition each sample on it** — arm A of the step-9 seam, and
/// the whole of `doc/devel/ng/spec/calling_em_loop.md` §2's pseudocode in one function.
///
/// Three loops, one inside the next, and **ng ships with the outer two switched off**:
///
/// 1. the **discovery round** of spec §2's pseudocode (§4.1) — structurally a loop here, and
///    permanently a single pass: milestone E1 found the eligible set reads no posterior, so
///    discovery's built setting runs as a **pre-pass inside candidate selection** before this
///    arm is ever called (`DiscoveryMode::BeforeTheLoop`;
///    `doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §6.5), and the loop it
///    reaches is this one, over a table already widened;
/// 2. the **slippage round**, which re-fits this locus's slippage numbers from its own
///    reads and rebuilds the table under them (§5.1) — built, and `max_rounds = 0` by
///    default, which is the frozen setting: the round's body never runs;
/// 3. the **frequency loop**, the innermost, and the only one that repeats at the shipped
///    configuration ([`run_frequency_loop`]).
///
/// **Discovery's loop is written with a body that runs once**: its body ends in an
/// unconditional `break`, and that — not the configuration — is what makes it run once. E1's
/// pre-pass finding made that permanent rather than provisional — a round here would look at
/// the same evidence selection already looked at and admit nothing, which an E1 test asserts.
/// Validation is the second lock: the two round-wrapped discovery modes are refused before
/// they can reach a [`RunnableCallingLoopConfig`], so no run arrives here expecting rounds it
/// will not get. **The slippage round is the one with a live body**: it reads
/// `config.slippage_refit`, and at zero rounds — the shipped default — its body is one
/// comparison and a `break`, with no table rebuilt and no read attributed.
///
/// **The bodies are not empty — they hold the whole of the work.** Computing the locus's
/// emissions, assembling the likelihood table from them and running the frequency loop are
/// inside the inner one; what a round adds is a reason to go round again, not the work.
///
/// **On the SNP/indel path both rounds are ignored structurally rather than half-honoured**
/// (spec §5.1's closing paragraph): there are no slippage numbers to re-fit, and discovery's
/// retrace is defined on stutter attribution. A run that switched the re-fit on therefore
/// re-fits its tracts and leaves its ordinary sites untouched, with no per-locus setting
/// anywhere.
///
/// **At the defaults the locus's emissions are therefore computed exactly once**, before the
/// frequency loop, whatever the pass count — which is what makes the expensive half of the
/// likelihood's cost independent of it (§2, §8), and is what
/// [`EmissionCost`](crate::ng::calling::EmissionCost) records. **The genotype-likelihood table
/// assembled from them is written once too where nothing is contaminated**, because with no
/// fraction fitted the assembled row reads no allele frequency; where something is, it is
/// written again at the head of every pass and once more against the frequencies the loop
/// settled on (`doc/devel/ng/spec/read_likelihoods.md` §3.6).
///
/// # A sample the candidate step ruled uncallable leaves before the first pass
///
/// It is given **no scratch row at all** (spec §5.0's ruling, and the choice B2 and C3b both
/// recorded as this step's). The rows are the run's sample order with the gaps closed up, so
/// the M-step sums the rows there are, the convergence delta divides by the chromosomes those
/// rows carry, and the site quality's count axis runs over the same cohort — none of them
/// having to be told to skip anything. The final pass walks the run's samples against a row
/// cursor and writes [`SampleGenotypeCall::Missing`] where a sample has no row.
///
/// # The type parameter, and why a SNP/indel locus carries it
///
/// It is the repeat-tract emission model — the seam a bake-off between two of them needs. A
/// SNP/indel locus never reaches it, and the arm still names it, because one worker's scratch
/// is typed by it and the same arm calls both paths.
pub struct SummariseConditionLoop<Model, Prior> {
    emission: Model,
    prior: Prior,
}

impl<Model, Prior> SummariseConditionLoop<Model, Prior> {
    /// Arm A, scoring repeat tracts with `emission` and every locus's genotypes under `prior`.
    ///
    /// **Both are values rather than constants**, and for the same reason: each is a seam the
    /// design exists to compare across. The genotype prior's two implementations disagree by
    /// 11 points of genotype accuracy on GIAB at 5× (`spec/calling_priors.md` §2.2), and a run
    /// has to be able to say which it used.
    pub fn new(emission: Model, prior: Prior) -> Self {
        Self { emission, prior }
    }
}

impl<Model, Prior> LocusGenotyper<Model::Scratch> for SummariseConditionLoop<Model, Prior>
where
    Model: SsrEmissionModel,
    Prior: GenotypePriorModel,
{
    fn name(&self) -> &'static str {
        "summarise the cohort, then condition each sample on it (arm A)"
    }

    #[expect(
        clippy::never_loop,
        reason = "the discovery round of spec §2's pseudocode is structurally present and \
                  permanently a single pass: discovery's built setting is a pre-pass in \
                  candidate selection (E1's finding — the eligible set reads no posterior), \
                  so no mode fills this body, and its unconditional `break` is what makes it \
                  run once. The slippage round inside it genuinely loops now — its body is \
                  built — so this expectation is the outer loop's alone"
    )]
    fn call_locus(
        &self,
        evidence: &LocusEvidence<'_>,
        parameters: &FrozenParameters<'_>,
        candidates: CandidateAlleles,
        config: &RunnableCallingLoopConfig,
        scratch: &mut CallingScratch<Model::Scratch>,
    ) -> LocusInference {
        evidence.assert_matches_locus_and_run(&candidates, parameters);
        // **This arm turns no locus away at its front door, and it used to turn away two
        // classes.** A contaminated run was refused here until the SNP/indel mixture was built,
        // and every repeat tract until the tract branches below were; a contaminated *tract*
        // was refused one module away until the third term of its mixture was. All three are
        // now called. What remains are checks on whether the caller's tables describe one
        // locus, and those belong where the tables meet rather than here.
        let table = GenotypeTable::build(parameters.ploidy(), candidates.len());
        let genotypes = table.view();

        let callable_sample_count = evidence.callable_sample_count();
        assert!(
            callable_sample_count > 0,
            "every one of the {} samples at {} was ruled uncallable by the candidate step, so \
             this locus has nobody to call. **This is a precondition of this trait**, not a \
             statement about the data: a caller asks `LocusEvidence::callable_sample_count` \
             first and counts such a locus rather than offering it here (owner's ruling, \
             2026-09-01), because the scratch below cannot be prepared for no rows",
            evidence.sample_count(),
            evidence.region()
        );
        // **Prepared first, claimed second, and the order is checked rather than merely
        // conventional.** `prepare_for_locus` clears the row map, so claiming before it would
        // throw the claims away and leave the locus with none — which the first read of the
        // map then refuses by name, "0 of them were claimed" (measured by swapping the two).
        scratch.prepare_for_locus(callable_sample_count, &candidates, &genotypes);
        for run_sample in 0..evidence.sample_count() {
            if is_callable(evidence, run_sample) {
                scratch.claim_row_for(
                    run_sample,
                    parameters.inbreeding_coefficient_by_sample()[run_sample],
                );
            }
        }
        // **Whether this locus's genotype-likelihood table moves as the loop iterates**, which
        // is not the same question as whether the run is contaminated.
        //
        // At an **ordinary site** it is: `q(o)`, the contaminating population's frequency for
        // the allele an observation shows, is the cohort's own estimate at this locus and the
        // loop rewrites it at every pass (`spec/read_likelihoods.md` §3.6). At a **repeat
        // tract** it is not: the third term there is the *fit's* length spectrum for this
        // tract's stratum, frozen before calling starts — §4.5.1 weighed the cohort's own
        // frequencies against it and refused them, because contamination must not move from one
        // pass to the next. So a contaminated tract is scored from one table like any other.
        let contaminated = !parameters.contamination_is_absent();
        let table_moves_with_the_loop =
            contaminated && matches!(evidence, LocusEvidence::Generic { .. });
        // **The contaminant tables are sized only where there is a mixture to fill them**, and
        // `prepare_for_locus` un-sized them a moment ago — so an uncontaminated locus cannot
        // read a contaminated one's frequencies, and a contaminated one cannot be scored
        // against tables nobody prepared. A tract needs none of them: nothing it scores reads a
        // per-batch frequency.
        if table_moves_with_the_loop {
            scratch.prepare_contaminant_tables(parameters.batch_count(), parameters.sample_count());
        }

        // **The locus's seed concentration: what the prior behaves as though it had already
        // seen here, before any sample's reads** (`spec/calling_priors.md` §2.3). Once per
        // locus — no pass moves it. Each fill returns a view over the buffer it just wrote;
        // what every pass reads is the buffer, through the scratch, so the view is dropped
        // here.
        //
        // **The two paths seed from two different measurements, and the specs keep the words
        // apart** (`spec/population_diversity.md` §2). An ordinary site is seeded from the
        // population's **frequency spectrum** — how allele frequencies are spread across the
        // population — which the run projected once, before any locus, into the two numbers and
        // the regime mark a `SpectrumSeed` carries. A repeat tract is seeded from its stratum's
        // **length spectrum** — how that stratum's chromosomes are spread over tract lengths —
        // which is looked up per locus.
        //
        // **Both have a fallback ladder and only one of them reports which rung it used.** The
        // ordinary-site ladder is `population_diversity.md` §3.4's, as amended by
        // `ordinary_site_seed.md` §4, and its rung is resolved when the run's seed is projected
        // and is not carried onto a locus; the tract ladder's three rungs are §4.4's, and the
        // rung comes back on the lookup and travels onto this locus's record.
        let length_spectrum_rung = match evidence {
            LocusEvidence::Generic { .. } => {
                let _ = fill_locus_concentration(
                    parameters.prior_seed(),
                    variant_class_of(&candidates),
                    candidates.len(),
                    scratch.seed_concentration_mut(),
                );
                None
            }
            LocusEvidence::Ssr {
                detail,
                candidate_repeat_counts,
                ..
            } => {
                // **The reference tract's own repeat count is entry 0**, and that is a
                // property of the candidate table rather than a convention agreed here: the
                // reference allele is id 0 of every `CandidateAlleles` and stays there, and
                // these counts are parallel to it. It is what the length spectrum's offsets are
                // measured from, and passing a *candidate's* count instead would re-centre
                // the shape on that candidate and flatten the prior.
                //
                // **A run whose fit produced no length spectrum anywhere is called here, at
                // the ladder's bottom rung, rather than refused** — the owner's ruling of
                // 2026-08-27, taken where `population_diversity.md` §5's *refuse a tract in a
                // run carrying no tract parameters* met §4.4's *the ladder always answers*.
                // **Refusing would turn a whole class of runs into a hard failure for a
                // condition the output already states**: the rung on the locus's record says
                // `StatedFlat`, which is how a call resting on a stated constant stays
                // distinguishable from one resting on a measurement without re-running
                // anything. §5 and §6 of that spec carry the ruling.
                let reference_repeats = candidate_repeat_counts[0].get();
                let length_spectrum = parameters.ssr_length_spectrum_at(
                    detail.motif.ssr_period(),
                    RepeatCount(reference_repeats),
                );
                // **One vector of `candidates` counts per tract**, and the third the repeat
                // path allocates — `fill_tract_likelihood_table` names the other two. It is a
                // type conversion: the seed builder takes plain counts, and the evidence
                // carries them as non-zero ones, which is the stronger type and the one worth
                // keeping at the seam.
                let repeats_of_each_candidate: Vec<u32> = candidate_repeat_counts
                    .iter()
                    .map(|count| count.get())
                    .collect();
                let _ = fill_ssr_seed(
                    &repeats_of_each_candidate,
                    reference_repeats,
                    length_spectrum,
                    scratch.seed_concentration_mut(),
                );
                Some(length_spectrum.rung())
            }
        };

        // **Whether this locus runs the slippage round at all** (spec §5.1). Frozen —
        // `max_rounds` 0, the shipped configuration — is this code at zero rounds: the test
        // below is the only thing that executes for it, no table is rebuilt and no read is
        // attributed. An ordinary site ignores the setting structurally, whatever a run asked
        // for: it has no slippage numbers to re-fit (§5.1's closing paragraph).
        let refit_config = &config.slippage_refit;
        let refit_may_run =
            !refit_config.is_frozen() && matches!(evidence, LocusEvidence::Ssr { .. });
        // The round's working memory, allocated only where a round may run: the emissions the
        // build caches for the attribution, the numbers the last round adopted, and a buffer
        // for the next re-fit. The frozen path allocates none of it — the three empty
        // constructors below hold no heap until something pushes into them.
        let mut cached_emissions = refit_may_run.then(RefitEmissionCache::default);
        let mut adopted_cells: Vec<Option<Slippage>> = Vec::new();
        let mut refitted_cells: Vec<Option<Slippage>> = Vec::new();
        let mut rounds_adopted: u32 = 0;

        // ── the discovery round (spec §4.1), `Off` by default ──
        let outcome = loop {
            // ── the slippage round (spec §5.1), `max_rounds` 0 by default ──
            let outcome = loop {
                // **The expensive half, once per round**: the emissions read no allele
                // frequency, so nothing the frequency loop does can move them — only a
                // slippage round that adopted new numbers can, by coming back around here.
                fill_what_no_pass_recomputes_for_round(
                    &self.emission,
                    evidence,
                    parameters,
                    &candidates,
                    &genotypes,
                    scratch,
                    SlippageRoundBuild {
                        adopted_slippage: (rounds_adopted > 0).then_some(adopted_cells.as_slice()),
                        emissions: cached_emissions.as_mut(),
                    },
                );
                // **And the cheap half, once here and then once per pass where something is
                // contaminated.** This first one is what the prior-free initialisation pass
                // reads, and there is no estimate of the frequencies for it to score against —
                // no pass has run — so it scores the reads alone, which is that pass's whole
                // job (§3).
                let reassembly =
                    table_moves_with_the_loop.then(|| TableReassembly::of(evidence, parameters));
                assemble_genotype_likelihood_table(
                    evidence,
                    parameters,
                    &genotypes,
                    if contaminated {
                        ContaminantFrequencies::TheReadsAlone
                    } else {
                        ContaminantFrequencies::NoneFitted
                    },
                    scratch,
                );
                let outcome = run_frequency_loop(
                    scratch,
                    &genotypes,
                    &self.prior,
                    config,
                    reassembly.as_ref(),
                );
                // **One more assembly against the frequencies the loop settled on**, so that
                // the final pass scores every sample against the same estimate the convergence
                // test just accepted rather than against the one the last pass started from.
                // Where nothing is contaminated the table reads no frequency and there is
                // nothing to assemble again.
                //
                // **What it costs is worth naming**: the genotypes and the site quality then
                // come from a table one assembly *newer* than the one the convergence test
                // looked at. At a locus that settled, the two differ by less than the
                // convergence threshold by construction; at one that hit the pass cap, the
                // difference is whatever the last pass moved — and such a locus is emitted
                // flagged, which is what that flag is for (§6).
                if let Some(reassembly) = reassembly.as_ref() {
                    reassembly.assemble(
                        &genotypes,
                        ContaminantFrequencies::TheLoopsOwnEstimate,
                        scratch,
                    );
                }
                // ── the round's own body: attribute, re-fit, and either stop or go round ──
                //
                // The frozen arm and every ordinary site break here having built one table,
                // which is the whole of the frozen behaviour.
                if !refit_may_run {
                    break outcome;
                }
                // The cap counts **adopted** rounds — each one a table rebuild and a rerun of
                // the frequency loop — mirroring production's `for _ in 0..refit_max_rounds`
                // (`src/ssr/cohort/em.rs:572`): at the cap, the numbers the last round
                // adopted stand and no further attribution runs.
                if rounds_adopted >= refit_config.max_rounds {
                    break outcome;
                }
                let tract = tract_evidence_of(evidence);
                let pooled = attribute_reads_under_the_posteriors(
                    scratch,
                    &genotypes,
                    &self.prior,
                    &tract,
                    &candidates,
                    cached_emissions
                        .as_ref()
                        .expect("a round that may run cached its build's emissions"),
                );
                refit_cells(
                    scratch.tract_fits().frozen_slippage_cells(),
                    &pooled,
                    refit_config,
                    &mut refitted_cells,
                );
                // **Test before adopt, and a converging round's numbers are discarded** —
                // production's shape exactly (`em.rs:576-580`): the round that moves less
                // than the threshold keeps the numbers the table was already built under,
                // so the emitted calls were scored under the numbers the test accepted.
                let standing: &[Option<Slippage>] = if rounds_adopted == 0 {
                    scratch.tract_fits().frozen_slippage_cells()
                } else {
                    &adopted_cells
                };
                if largest_movement(standing, &refitted_cells)
                    < refit_config.round_convergence_threshold
                {
                    break outcome;
                }
                std::mem::swap(&mut adopted_cells, &mut refitted_cells);
                rounds_adopted += 1;
                // And round again: rebuild the table under the adopted numbers and rerun
                // the frequency loop to convergence — the nested shape spec §5.1 decides.
            };
            // Discovery off: the candidate set is what selection settled on, and the prune
            // and the re-run that would follow a round have nothing to do.
            break outcome;
        };
        // The probe for the round's stopping rule — nothing on the scoring path reads it.
        scratch.record_slippage_refit_rounds(rounds_adopted);

        // **Read before the final pass rather than inside its argument list**, because at a
        // tract it reads the scratch's own parameter table and the pass takes that scratch
        // mutably.
        let weakest_warrant =
            weakest_warrant_at_the_locus(evidence, parameters, scratch.tract_fits());
        // **What a repeat tract's call rested on, read off the fits the gather just filled** —
        // `None` at an ordinary site, which is what keeps the two paths' ladders apart:
        // `LocusInference::new` refuses a record at a SNP or indel and refuses its absence at a
        // tract.
        //
        // **Built from the rung rather than from the fits' own state**, and that is what makes
        // it safe to read a per-worker buffer here: the fits are cleared and refilled per tract
        // and hold the *last* tract's counts at an ordinary site, so asking them whether they
        // describe this locus would answer about a locus already emitted. The rung is `Some`
        // exactly at a tract and was resolved for this one.
        let repeat_tract = length_spectrum_rung.map(|rung| {
            let fits = scratch.tract_fits();
            RepeatTractProvenance::new(
                rung,
                fits.cell_count(),
                fits.cells_with_no_fitted_slippage(),
                fits.cells_whose_read_group_the_fit_does_not_describe(),
                fits.cells_with_no_fitted_substitution_rate(),
                fits.contaminant_term_was_built(),
            )
        });
        summarise_final_pass(
            scratch,
            &genotypes,
            evidence,
            parameters,
            &self.prior,
            candidates,
            outcome,
            weakest_warrant,
            repeat_tract,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::parameter_estimation::joint::stratum_fits::LengthSpectrumRung;

    /// The repeat unit every repeat tract in this file is built on — **a dinucleotide**, so that
    /// a whole repeat and a single base are different steps and the part-repeat half of the
    /// stutter distribution is not hidden.
    const TRACT_MOTIF: &[u8] = b"AT";

    /// The two candidate lengths every repeat tract here is called over, in candidate-table
    /// order: **6 whole repeats and 7**. Entry 0 is the reference allele's, which is what the
    /// prior's length spectrum measures its offsets from.
    ///
    /// **Different from each other**, so a lookup that read one candidate's repeat count for
    /// another lands in a different stratum rather than the same one. **Neither is what a
    /// reader could derive from the other's bases**, which
    /// `an_interrupted_candidates_supplied_repeat_count_is_what_picks_its_stratum` is about.
    const TRACT_CANDIDATE_REPEATS: [u32; 2] = [6, 7];

    /// **How many libraries a tract fixture's run has, and it is not two.** The parameter table
    /// a tract is scored from is `read groups × candidates`, filled read-group-major; at an
    /// equal shape a table filled the other way round is the same length and the same set of
    /// cells, so a transposition passes every shape check. At three libraries against two
    /// candidates the two orders disagree cell by cell.
    ///
    /// Only the first library sends any read below. The other two are in the run, which is what
    /// makes them visible to a tract's warrant fold — that fold runs over the run's read-group
    /// axis rather than over the reads that arrived.
    const TRACT_READ_GROUPS: usize = 3;

    /// The two candidates' repeat counts. **Supplied, not selected** — repeat-tract candidate
    /// selection is unwritten, so a tract fixture states its candidates' repeat counts and a
    /// reader must not take them for a step's output.
    fn tract_repeat_counts() -> Vec<NonZeroU32> {
        TRACT_CANDIDATE_REPEATS
            .iter()
            .map(|count| NonZeroU32::new(*count).expect("a candidate always holds a repeat"))
            .collect()
    }

    /// The bases of a tract holding `repeats` whole copies of the motif.
    fn tract_bases(repeats: u32) -> Vec<u8> {
        TRACT_MOTIF.repeat(repeats as usize)
    }

    fn tract_detail() -> SsrDetail {
        SsrDetail {
            motif: Motif::new(TRACT_MOTIF).expect("a dinucleotide motif"),
            left_flank: Box::from(b"CCCGGG".as_slice()),
            right_flank: Box::from(b"TTTAAA".as_slice()),
        }
    }

    /// The candidate table the two lengths above spell, with the shorter as the reference.
    fn tract_alleles() -> CandidateAlleles {
        let mut alleles = CandidateAlleles::new(
            tract_bases(TRACT_CANDIDATE_REPEATS[0]).into_boxed_slice(),
            LocusKind::Ssr(tract_detail()),
        );
        alleles.admit(tract_bases(TRACT_CANDIDATE_REPEATS[1]).into_boxed_slice());
        alleles
    }

    /// `reads` reads that spanned the whole tract and showed `repeats` whole copies of it, from
    /// the run's first library.
    ///
    /// **The per-read error sum is set to a plausible value and nothing on this path reads
    /// it**, which is worth saying because it is what a reader would assume decides these
    /// calls. A tract's row is scored from a stutter model and a per-base substitution rate;
    /// `q_sum` is the SNP/indel emission's charge and neither `likelihood::ssr` nor its
    /// emission model touches it. Zeroing it here would change no number below.
    fn tract_reads(repeats: u32, reads: u32) -> SequenceObservation {
        SequenceObservation {
            bases: tract_bases(repeats).into_boxed_slice(),
            read_witness: ReadWitness::Complete,
            read_group: ReadGroupId(0),
            num_obs: reads,
            num_fwd: reads / 2,
            q_sum: SummedLogError::from_nats(-10.0 * f64::from(reads)),
            mapq_sum: 60 * reads,
            mapq_sum_sq: u64::from(reads) * 3_600,
            placed_left: reads / 2,
            chain_ids: Vec::new(),
        }
    }

    /// **`reads` reads that ran out inside the tract**, having seen `repeats` whole copies of it
    /// before they did — a lower bound on the tract's length rather than a measurement of it.
    ///
    /// **Over half the reads overlapping a 60-base tract are these on real data**, and they are
    /// scored by a different term of the row from a read that spanned the whole of it. A fixture
    /// of complete reads only is the one shape in which counting the observations and counting
    /// the *complete* observations are the same number.
    fn tract_reads_that_ran_out(repeats: u32, reads: u32) -> SequenceObservation {
        let bases = tract_bases(repeats);
        let covered = u16::try_from(bases.len()).expect("a short fixture tract");
        SequenceObservation {
            read_witness: ReadWitness::from_left(
                covered,
                crate::ng::locus_generation::LocusLen::from_positions(u64::from(covered) + 4),
            )
            .expect("a reach that covers part of the locus"),
            bases: bases.into_boxed_slice(),
            read_group: ReadGroupId(0),
            num_obs: reads,
            num_fwd: reads / 2,
            q_sum: SummedLogError::from_nats(-10.0 * f64::from(reads)),
            mapq_sum: 60 * reads,
            mapq_sum_sq: u64::from(reads) * 3_600,
            placed_left: reads / 2,
            chain_ids: Vec::new(),
        }
    }

    /// **A repeat fit that reached every stratum this file's tracts are called over**, each
    /// with its own slippage numbers and its own length spectrum.
    ///
    /// **Three sources of coincidence are removed on purpose**, each because a mutation
    /// survived it:
    ///
    /// - **the two strata share no number.** Level, direction split and fall-off all differ, so
    ///   a lookup that answered from the wrong stratum changes a score rather than only the one
    ///   number a reader happens to check;
    /// - **the run's libraries do not share one slippage group.** Library 0 sits in group 0 and
    ///   libraries 1 and 2 in group 1, whose numbers differ — so a parameter table filled
    ///   candidate-major instead of read-group-major hands library 0's reads another group's
    ///   polymerase;
    /// - **no length spectrum is a palindrome, no two of its classes share a weight, and no two
    ///   adjacent pairs of classes are in the same ratio.** Each is heavier on the contraction
    ///   side, as every real repeat fit is, so reversing one is a different prior; reading class
    ///   *i* for class *j* is a different number; and — the repair a review forced — a spectrum
    ///   read one repeat off centre gives the candidates a *different* pair of shares, where a
    ///   geometric tail would have given the same pair.
    ///
    /// Every warrant here is the cell's own, so a tract scored under it is
    /// [`Provenance::FittedHere`] — which is what makes the calibration's own `Defaulted`
    /// visible if it ever entered the tract's fold.
    fn tract_strata() -> StratumFits {
        tract_strata_describing(std::collections::BTreeMap::from([
            (ReadGroupId(0), 0),
            (ReadGroupId(1), 1),
            (ReadGroupId(2), 1),
        ]))
    }

    /// The same fit, over whichever of the run's libraries it claims to describe — what the
    /// fixture about a library the fit never saw needs.
    fn tract_strata_describing(
        slippage_group_of_each_library: std::collections::BTreeMap<ReadGroupId, u32>,
    ) -> StratumFits {
        use crate::ng::parameter_estimation::joint::census::Stratum;
        use crate::ng::parameter_estimation::joint::share_curve::ShareSource;
        use crate::ng::parameter_estimation::joint::slippage_curve::LevelSource;
        use crate::ng::parameter_estimation::joint::ssr_fit::{
            LevelProvenance, ShareProvenance, SharesProvenance, Slippage, StratumFit,
            StratumOutcome,
        };

        let level = LevelProvenance {
            source: LevelSource::Cell,
            curve: None,
            reach: None,
            slipped_reads: Some(400.0),
        };
        let share = ShareProvenance {
            source: ShareSource::Stratum,
            curve: None,
            reach: None,
        };
        // One entry per slippage group: group 0 is library 0's, group 1 is libraries 1 and 2's.
        let slippage_of_each_group = |level_in_group_0: f64| {
            vec![
                Some(Slippage {
                    level: level_in_group_0,
                    shorter_share: 0.83,
                    fall_off: 0.25,
                }),
                Some(Slippage {
                    level: level_in_group_0 * 1.75,
                    shorter_share: 0.62,
                    fall_off: 0.41,
                }),
            ]
        };
        let fitted_stratum =
            |repeats: u64, level_in_group_0: f64, length_spectrum: Vec<f64>, concentration: f64| {
                StratumOutcome::Fitted(Box::new(StratumFit {
                    stratum: Stratum {
                        period: 2,
                        reference_repeats: repeats,
                    },
                    slippage: slippage_of_each_group(level_in_group_0),
                    length_spectrum,
                    concentration,
                    log_likelihood_a_tract: -1.5,
                    tracts_fitted: 40,
                    borrowed: Vec::new(),
                    converged: true,
                    tracts_of_its_own: 40,
                    reads_crossing: 400,
                    level_provenance: vec![Some(level), Some(level)],
                    shares_provenance: vec![
                        Some(SharesProvenance {
                            slipped_reads: Some(400.0),
                            shorter_share: share,
                            fall_off: share,
                        }),
                        Some(SharesProvenance {
                            slipped_reads: Some(400.0),
                            shorter_share: share,
                            fall_off: share,
                        }),
                    ],
                }))
            };
        StratumFits::over(
            &[
                // The weights run over offsets −2 … +2 whole repeats from the stratum's own
                // reference tract length.
                //
                // **No two adjacent pairs of these weights are in the same ratio**, which is a
                // deliberate repair: the first version ran `0.10, 0.25, 0.45, 0.15, 0.05`,
                // whose upper tail falls by a factor of three at each step — so a spectrum read
                // one repeat off centre gave the two candidates the *same* pair of shares,
                // 0.75 and 0.25, and the re-centring this fixture exists to catch was
                // invisible. Here offsets {0, +1} give 0.8 and 0.2 where {+1, +2} give 0.6875
                // and 0.3125.
                fitted_stratum(4, 0.03, vec![0.12, 0.27, 0.42, 0.14, 0.05], 16.0),
                fitted_stratum(6, 0.04, vec![0.10, 0.30, 0.44, 0.11, 0.05], 20.0),
                fitted_stratum(7, 0.06, vec![0.09, 0.21, 0.43, 0.19, 0.08], 25.0),
            ],
            slippage_group_of_each_library,
        )
    }

    /// **The record a called repeat tract carries**, which every tract has and no ordinary site
    /// does — so a call to this on a SNP/indel locus's result is a test asking the wrong path.
    fn tract_record(inference: &LocusInference) -> RepeatTractProvenance {
        inference
            .repeat_tract
            .expect("every called repeat tract carries what its parameters rested on")
    }

    /// Which rung of the tract ladder that record's prior shape came from.
    fn rung_of(inference: &LocusInference) -> LengthSpectrumRung {
        tract_record(inference).length_spectrum_rung()
    }

    /// The run's libraries, none of them calibrated — a tract's row never reads a calibration,
    /// so this says nothing about the tract and everything about whether the warrant fold has
    /// wandered onto the wrong list.
    fn tract_libraries() -> Vec<ReadGroupCalibration> {
        vec![ReadGroupCalibration::defaulted(); TRACT_READ_GROUPS]
    }

    /// **A fitted per-base substitution rate for every `(library, stratum)` cell a tract's
    /// parameter table covers** — every library of the run against every stratum this file's
    /// tracts reach.
    ///
    /// **Every cell, not only the ones whose reads arrived**, because the table is built on the
    /// run's read-group axis: a run that fitted only the library that sent reads still has
    /// cells for the other two, and those cells fall to the stated constant and pull the
    /// locus's warrant down with them. `a_library_that_sent_this_tract_no_read_still_weakens_it`
    /// is that case; this is its opposite.
    ///
    /// The rate differs per stratum, so a lookup keyed by the wrong repeat count is a different
    /// number rather than the same one.
    fn tract_substitution_rates() -> std::collections::BTreeMap<
        crate::ng::parameter_estimation::ssr::StratumKey,
        crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
    > {
        substitution_rates_over(&[4, 6, 7], TRACT_READ_GROUPS)
    }

    /// The same, over the strata and libraries named — so a fixture can leave a cell out.
    fn substitution_rates_over(
        strata: &[u32],
        read_groups: usize,
    ) -> std::collections::BTreeMap<
        crate::ng::parameter_estimation::ssr::StratumKey,
        crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
    > {
        use crate::ng::parameter_estimation::Estimate;
        use crate::ng::parameter_estimation::ssr::{
            RepeatCount, Stratum as SsrStratum, StratumKey,
        };
        use crate::ng::types::{ErrorRate, SsrPeriod};

        let period = SsrPeriod::try_new(2).expect("a dinucleotide");
        let mut rates = std::collections::BTreeMap::new();
        for library in 0..read_groups {
            for repeats in strata {
                rates.insert(
                    StratumKey {
                        read_group: ReadGroupId(
                            u32::try_from(library).expect("a small library index"),
                        ),
                        stratum: SsrStratum::new(period, RepeatCount(*repeats)),
                        ploidy: diploid(),
                    },
                    Estimate {
                        // A different rate per stratum, so a lookup keyed by the wrong repeat
                        // count is a different number.
                        value: ErrorRate::try_new(0.001 * f64::from(*repeats))
                            .expect("a probability"),
                        provenance: Provenance::FittedHere,
                        observations: 4_000,
                    },
                );
            }
        }
        rates
    }
    use proptest::prelude::*;

    use crate::ng::calling::GenericObservation;
    use crate::ng::calling::genotype_prior::{
        MarginalizedDirichletPrior, SeedRegime, SpectrumSeed,
    };
    use crate::ng::calling::likelihood::ssr_emission::SsrCandidate;
    use crate::ng::calling::quality::MAX_GENOTYPE_QUALITY;
    use crate::ng::calling::{
        CandidateAlleles, ExpectedAlleleCopies, FrozenParameters, GenotypeIdx, GenotypeTable,
        LocusEvidence, LocusInference, ReadGroupCalibration, SampleGenotypeCall, SsrSampleEvidence,
        UNWRITTEN_SCRATCH_VALUE,
    };
    use crate::ng::locus_generation::{
        LocusKind, ReadWitness, SequenceObservation, SsrDetail, WitnessedLocusPositions,
    };
    use crate::ng::parameter_estimation::Provenance;
    use crate::ng::parameter_estimation::joint::stratum_fits::StratumFits;
    use crate::ng::run::cohort_merge::build::PartialObservation;
    use crate::ng::types::{
        AlleleId, ContigId, GenomeRegion, Genotype, Motif, Phred, Position, ReadGroupId,
        SummedLogError,
    };
    use std::num::NonZeroU32;
    use std::sync::Arc;

    use super::super::{CallingLoopConfig, DEFAULT_CONVERGENCE_THRESHOLD};
    use crate::ng::calling::MIN_CONTAMINANT_FREQUENCY;
    use crate::ng::calling::tests::one_batch;
    use crate::ng::parameter_estimation::joint::sequencing_batches::SequencingBatches;
    use crate::ng::read::input::read_groups::ReadGroups;
    use crate::ng::types::BatchId;

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

    /// A run whose repeat-tract substitution rates were never fitted — the empty map, which is
    /// what a fixture calling no tract needs, and what
    /// `FrozenParameters::ssr_substitution_rate_at` answers `None` from.
    ///
    /// A `static` rather than a function, so that a call site can borrow it for as long as the
    /// parameters live: `BTreeMap::new` is a `const fn`, and a temporary would be freed at the
    /// end of the statement that built the view.
    static NO_SUBSTITUTION_RATES: std::collections::BTreeMap<
        crate::ng::parameter_estimation::ssr::StratumKey,
        crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
    > = std::collections::BTreeMap::new();

    /// A site quality standing in for one the worker computed, where what the test is about
    /// is something else.
    ///
    /// **Deliberately not zero.** Zero is a real answer — a locus at which nobody carries
    /// anything — so a test whose fixture wrote zero would pass against an implementation
    /// that lost the value on the way in.
    fn a_worker_written_site_quality() -> Phred {
        Phred::try_new(37.0).expect("a legal quality, and below the site-quality ceiling")
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
    /// winners off pass 1 passes against an M-step that computes nothing. Measured during
    /// B2's review: with the summing loop deleted, the one-pass version of the permutation
    /// test below stayed green while three other tests failed. From pass 2 the winners are
    /// scored against the row the previous pass's M-step produced, which is the regime the
    /// coupling under test actually occurs in.
    ///
    /// The posterior is one reused buffer, so each sample's winner is read off before the next
    /// is scored — which is also why step C3 of the plan has to take the genotype quality as it
    /// scores rather than afterwards.
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
    /// becomes permanent is governed by the pass cap, which ships at 50 and is C2's, and by
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
    /// Before C1 this width was checked only inside `PriorRow::new`, which a flat pass never
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

    /// **A sample with no reads votes for a 50% allele frequency on the flat pass**, and the
    /// design does not currently say whether it should.
    ///
    /// `spec/calling_em_loop.md` §7 says a sample with no coverage *"scores every genotype
    /// alike, so the prior decides it alone — the right answer rather than a special case"*.
    /// **On a flat pass there is no prior to decide it**, so its posterior is the normalised
    /// genotype table and its expected copies come out as the average genotype: a full copy of
    /// the alternative allele at a biallelic locus, the same contribution as a confident
    /// heterozygote.
    ///
    /// This test pins the behaviour rather than endorsing it — the number is `1.0`, and at
    /// three reads a position roughly one sample in twenty is silent at any given position, so
    /// this is the ordinary case rather than a corner. Raised for the owner against §3 and §7.
    #[test]
    fn a_sample_with_no_reads_contributes_a_full_alternative_copy_to_the_flat_pass() {
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

        let copies = scratch.sample_expected_copies(0);
        assert!(
            (copies[1] - 1.0).abs() < 1e-12,
            "a silent sample's flat-pass copies of the alternative allele: {copies:?}"
        );
        assert!(
            (copies[0] - 1.0).abs() < 1e-12,
            "and of the reference: {copies:?}"
        );
    }

    // ---------------------------------------------------------------------------------
    // C2 — convergence, the cap, and the emitted flag
    // ---------------------------------------------------------------------------------

    /// A configuration that ships, with the pass cap moved.
    fn capped_at(max_passes: u32) -> RunnableCallingLoopConfig {
        CallingLoopConfig {
            max_passes: NonZeroU32::new(max_passes).expect("a cap of at least one pass"),
            ..CallingLoopConfig::DEFAULT
        }
        .validate()
        .expect("only the cap moved, and it is not a value validation refuses")
    }

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

    /// Run the frequency loop over a cohort whose genotype likelihoods are `likelihoods`, one
    /// row per sample, every sample outbred — and hand back the scratch so a test can read
    /// what the loop left in it: the cohort's expected copies, the previous pass's row, and
    /// the last sample scored.
    ///
    /// For a cohort whose samples differ in inbreeding, use [`loop_over_inbred`].
    fn loop_over(
        likelihoods: &[Vec<f64>],
        seed_concentration: &[f64],
        alleles: &CandidateAlleles,
        view: &GenotypeTableView<'_>,
        config: &RunnableCallingLoopConfig,
    ) -> (FrequencyLoopOutcome, CallingScratch<()>) {
        let outbred_cohort = vec![outbred(); likelihoods.len()];
        loop_over_inbred(
            likelihoods,
            &outbred_cohort,
            seed_concentration,
            alleles,
            view,
            config,
        )
    }

    /// [`loop_over`] with an inbreeding coefficient chosen per sample, in the run's sample
    /// order.
    fn loop_over_inbred(
        likelihoods: &[Vec<f64>],
        inbreeding_by_sample: &[InbreedingF],
        seed_concentration: &[f64],
        alleles: &CandidateAlleles,
        view: &GenotypeTableView<'_>,
        config: &RunnableCallingLoopConfig,
    ) -> (FrequencyLoopOutcome, CallingScratch<()>) {
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(likelihoods.len(), alleles, view);
        // Every sample of this fixture is callable, so the rows are the run's samples
        // one for one — which is what a locus with nothing set aside looks like.
        for (sample, &inbreeding) in inbreeding_by_sample.iter().enumerate() {
            scratch.claim_row_for(sample, inbreeding);
        }
        scratch
            .seed_concentration_mut()
            .copy_from_slice(seed_concentration);
        for (sample, row) in likelihoods.iter().enumerate() {
            for (slot, &value) in scratch
                .sample_genotype_likelihoods_mut(sample)
                .iter_mut()
                .zip(row)
            {
                *slot = LogProb(value);
            }
        }
        let outcome = run_frequency_loop(
            &mut scratch,
            view,
            &MarginalizedDirichletPrior,
            config,
            None,
        );
        (outcome, scratch)
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

        let moved =
            (as_given.cohort_expected_copies()[1] - traded.cohort_expected_copies()[1]).abs();
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
        let _ = cohort_expected_copies_have_settled(
            &[],
            &[],
            diploid(),
            1,
            DEFAULT_CONVERGENCE_THRESHOLD,
        );
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

    // ────────────────────────────────────────────────────────────────────────────────
    // C3b — the final pass: the calls, the site quality, and the artifact summary
    // ────────────────────────────────────────────────────────────────────────────────

    fn locus_region() -> GenomeRegion {
        GenomeRegion {
            contig: ContigId(3),
            start: Position(940),
            end: Position(940),
        }
    }

    /// One `(allele, read group)` row of a sample's evidence, with the two counts the
    /// artifact summary reads.
    ///
    /// **`q_sum` is filled with something that would be wrong if anything read it**, rather
    /// than with zero: the final pass never looks at the error mass, and a fixture of zeroes
    /// would not say so.
    fn observation(
        allele: u16,
        read_group: u32,
        num_reads: u32,
        forward_reads: u32,
        placed_left_reads: u32,
    ) -> GenericObservation {
        GenericObservation {
            allele: AlleleId(allele),
            read_group: ReadGroupId(read_group),
            num_reads,
            q_sum: -3.0 * f64::from(num_reads),
            forward_reads,
            placed_left_reads,
        }
    }

    /// A sample whose reads are `rows` and which the candidate step did not set aside.
    fn called_sample<'a>(rows: &'a [GenericObservation]) -> GenericLocusSample<'a> {
        GenericLocusSample {
            evidence: GenericSampleEvidence::new(rows, 0.0, &[]),
            genotype_must_be_missing: false,
        }
    }

    /// A neutral panel's fitted spectrum at one variant per kilobase — human diversity, and
    /// the seed the site quality's own tests use.
    fn human_like_seed() -> SpectrumSeed {
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape)
    }

    /// **The whole of one locus over a hand-built likelihood table**: the frequency loop, then
    /// the final pass over the same scratch — the order and the sharing the design requires,
    /// since the pass reads the posterior row the loop leaves and the likelihood table beside
    /// it. **The table is supplied rather than assembled here**, so nothing in this fixture
    /// re-assembles it; the driver's own tests are where a contaminated locus's per-pass
    /// assembly is exercised.
    ///
    /// `likelihoods` holds **one row per sample of the run**, and a sample the candidate step
    /// ruled uncallable simply does not get a scratch row — so its row here is never read,
    /// which is what the driver does with it too.
    fn call_locus(
        likelihoods: &[Vec<f64>],
        evidence: &LocusEvidence<'_>,
        seed_concentration: &[f64],
        alleles: &CandidateAlleles,
        view: &GenotypeTableView<'_>,
        model: &dyn GenotypePriorModel,
    ) -> LocusInference {
        let calibration = [ReadGroupCalibration::defaulted()];
        let inbreeding = vec![outbred(); likelihoods.len()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let parameters = FrozenParameters::uncontaminated(
            &calibration,
            &inbreeding,
            human_like_seed(),
            &strata,
            &NO_SUBSTITUTION_RATES,
            view.ploidy(),
        );
        let callable: Vec<usize> = (0..evidence.sample_count())
            .filter(|run_sample| is_callable(evidence, *run_sample))
            .collect();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(callable.len(), alleles, view);
        for &run_sample in &callable {
            scratch.claim_row_for(run_sample, outbred());
        }
        scratch
            .seed_concentration_mut()
            .copy_from_slice(seed_concentration);
        for (row, &run_sample) in callable.iter().enumerate() {
            for (slot, &value) in scratch
                .sample_genotype_likelihoods_mut(row)
                .iter_mut()
                .zip(&likelihoods[run_sample])
            {
                *slot = LogProb(value);
            }
        }
        let outcome = run_frequency_loop(
            &mut scratch,
            view,
            model,
            &RunnableCallingLoopConfig::default(),
            None,
        );
        summarise_final_pass(
            &mut scratch,
            view,
            evidence,
            &parameters,
            model,
            alleles.clone(),
            outcome,
            Provenance::FittedHere,
            // **The fixture's own, because this helper drives the final pass and not the
            // loop**: nothing here gathers a tract's scoring parameters, so there is no filled
            // table to read the counts off. A tract must carry a record and an ordinary site
            // must not, so the locus's kind decides which.
            (!matches!(alleles.kind(), LocusKind::Generic))
                .then(a_tract_record_this_fixture_did_not_gather),
        )
    }

    /// **A stand-in for what a tract's parameters rested on**, for a fixture that never
    /// gathered any — the counts are the fixture's, not a measurement, and no test of this
    /// helper's may assert them.
    ///
    /// The driver's own tests are where the counts are checked against a gather.
    fn a_tract_record_this_fixture_did_not_gather() -> RepeatTractProvenance {
        RepeatTractProvenance::new(LengthSpectrumRung::StatedFlat, 0, 0, 0, 0, false)
    }

    /// The call this fixture expects, unwrapped — a `Missing` here is the test's failure and
    /// not something to match on.
    fn called(inference: &LocusInference, sample: usize) -> (&Genotype, Phred) {
        match &inference.per_sample[sample] {
            SampleGenotypeCall::Called {
                genotype,
                genotype_quality,
                ..
            } => (genotype, *genotype_quality),
            SampleGenotypeCall::Missing => {
                panic!("sample {sample} was called, so it is not missing")
            }
        }
    }

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
        let by_hand = 10.0 * (7.0_f64 / 3.0).log10();
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
    /// The two samples' calls do not enter the choice — it is made from the reads alone,
    /// before any sample is scored (`spec/calling_quality.md` §6.3).
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

    /// **Two alternatives the same number of reads reached: the lower id wins**, because the
    /// fold keeps the first strict maximum.
    ///
    /// The allele table's order is fixed, so this is the difference between a run that
    /// reproduces itself and one that does not — the same rule the genotype quality's argmax
    /// follows, for the same reason.
    #[test]
    fn a_tie_between_two_alternatives_goes_to_the_lower_allele_id() {
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

        let counts = inference
            .artifact_test_counts()
            .expect("both alternatives drew reads");
        assert_eq!(counts.primary_alternative, AlleleId(1));
        // And the counts are that allele's, not the other's — 3 forward against 1.
        assert_eq!(counts.alternative_forward_reads, 3.0);
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
            witnessed_in_locus: WitnessedLocusPositions::from_half_open_runs([(0_u16, 1_u16)])
                .expect("one witnessed position"),
            read_group: ReadGroupId(0),
            bases: Box::from(b"A".as_slice()),
            num_reads: 50,
            q_sum: -100.0,
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
        let per_sample = [SsrSampleEvidence::new(&[], &detail)];
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
        // compared against an independently prepared fold rather than merely being present.
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
            human_like_seed(),
        );
        assert!(expected.get() > 0.0, "the fixture's fold is above zero");
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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

    // ────────────────────────────────────────────────────────────────────────────────
    // D1 — the driver: the table built once, the outer rounds structurally off
    // ────────────────────────────────────────────────────────────────────────────────

    use crate::ng::calling::likelihood::ssr_emission::{
        StutterSubstitutionEmission, StutterSubstitutionScratch,
    };
    use crate::ng::calling::{ContaminationView, EmissionCost};
    use crate::ng::parameter_estimation::joint::contamination::ContaminationSource;

    /// The shipped arm: repeat tracts scored by the shipped emission model, genotypes under
    /// the marginalized Dirichlet prior.
    fn shipped_arm()
    -> SummariseConditionLoop<StutterSubstitutionEmission, MarginalizedDirichletPrior> {
        SummariseConditionLoop::new(StutterSubstitutionEmission, MarginalizedDirichletPrior)
    }

    /// A worker's scratch, typed by the shipped emission model.
    fn worker_scratch() -> CallingScratch<StutterSubstitutionScratch> {
        CallingScratch::default()
    }

    /// One run's frozen parameters over `samples` outbred samples and one read group, with
    /// nothing contaminated.
    fn uncontaminated_run<'a>(
        calibration: &'a [ReadGroupCalibration],
        inbreeding: &'a [InbreedingF],
        strata: &'a StratumFits,
        substitution: &'a std::collections::BTreeMap<
            crate::ng::parameter_estimation::ssr::StratumKey,
            crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
        >,
        ploidy: Ploidy,
    ) -> FrozenParameters<'a> {
        FrozenParameters::uncontaminated(
            calibration,
            inbreeding,
            human_like_seed(),
            strata,
            substitution,
            ploidy,
        )
    }

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
    /// builds, whatever the pass count, which is the property Milestone D exists for.
    #[test]
    fn the_driver_calls_genotypes_from_reads_and_builds_the_table_once() {
        let (alleles, _) = generic_locus(1);
        let carrier = [observation(1, 0, 8, 4, 4)];
        let reference_sample = [observation(0, 0, 8, 4, 4)];
        let per_sample = [called_sample(&carrier), called_sample(&reference_sample)];
        let evidence = LocusEvidence::generic(locus_region(), &per_sample);
        let calibration = [ReadGroupCalibration::defaulted()];
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
            (copies[0] - 2.015_067_414_196_071).abs() < 1e-12
                && (copies[1] - 1.984_932_585_803_929_3).abs() < 1e-12,
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
        let calibration = [ReadGroupCalibration::defaulted(); 3];
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
    /// every other test in the whole `ng::calling` module green, 733 of them, and moves this
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
        let calibration = [ReadGroupCalibration::defaulted(); 3];
        let fractions = [contaminated_at(0.05); 3];
        let inbreeding = [outbred(), outbred(), outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = SequencingBatches::all_together_over(3, 3);
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
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
        let copies = inference.cohort_expected_copies().copies();
        let expected = [3.008_738_983_972_081_7, 2.991_261_016_027_918];
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
            .zip([20.701_796_f32, 30.193_914, 17.712_873])
            .enumerate()
        {
            assert!(
                (got - want).abs() < 1e-4,
                "sample {sample}'s confidence: {got} against {want}"
            );
        }
        assert!(
            (inference.site_quality.get() - 119.721_42_f32).abs() < 1e-3,
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
        // confidence from **30.194 to 30.180** (measured 2026-08-26). The cohort's expected
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
            (capped.site_quality.get() - 119.720_436_f32).abs() < 1e-4,
            "the capped run's site quality: {}",
            capped.site_quality.get()
        );
        assert!(
            (called(&capped, 1).1.get() - 30.193_907_f32).abs() < 1e-4,
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
        let calibration = [ReadGroupCalibration::defaulted(); 3];
        let fractions = [contaminated_at(0.05); 3];
        let inbreeding = [outbred(), outbred(), outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = one_batch(3, 3);
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

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
        let calibration = [ReadGroupCalibration::defaulted(); 3];
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
        let calibration = [ReadGroupCalibration::defaulted(); 2];
        let fractions = [contaminated_at(0.05); 2];
        let inbreeding = [outbred(), outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = one_batch(2, 2);
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

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
                scratch.prepare_contaminant_tables(
                    parameters.batch_count(),
                    parameters.sample_count(),
                );
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
        let calibration = [ReadGroupCalibration::defaulted(); 4];
        let fractions = [contaminated_at(0.05); 4];
        let inbreeding = [outbred(), outbred(), outbred(), outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = two_batches_of_two();
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

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
        let calibration = [ReadGroupCalibration::defaulted(); 3];
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
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

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
        let calibration = [ReadGroupCalibration::defaulted(); 3];
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
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
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
        let calibration = [ReadGroupCalibration::defaulted()];
        let fractions = [contaminated_at(0.05)];
        let inbreeding = [outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = SequencingBatches::all_together_over(1, 1);
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

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
        let calibration = [ReadGroupCalibration::defaulted()];
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
            ReadGroupCalibration::defaulted(),
            ReadGroupCalibration::defaulted(),
        ];
        let fractions = [contaminated_at(0.05), contaminated_at(0.05)];
        let inbreeding = [outbred(), outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = one_batch_each();
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

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
            ReadGroupCalibration::defaulted(),
            ReadGroupCalibration::defaulted(),
        ];
        let fractions = [contaminated_at(0.05), contaminated_at(0.05)];
        let inbreeding = [outbred(), outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = one_batch_each();
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);

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
        let calibration = [ReadGroupCalibration::defaulted(); 2];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
            crate::ng::parameter_estimation::ssr::StratumKey,
            crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
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
            crate::ng::parameter_estimation::ssr::StratumKey,
            crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
        >,
        calibration: &[ReadGroupCalibration],
        scratch: &mut CallingScratch<StutterSubstitutionScratch>,
    ) -> LocusInference {
        let detail = tract_detail();
        let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
            .iter()
            .map(|observations| SsrSampleEvidence::new(observations, &detail))
            .collect();
        let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, repeat_counts);
        let inbreeding = vec![outbred(); per_sample.len()];
        let parameters =
            uncontaminated_run(calibration, &inbreeding, strata, substitution, diploid());
        shipped_arm().call_locus(
            &evidence,
            &parameters,
            alleles.clone(),
            &RunnableCallingLoopConfig::default(),
            scratch,
        )
    }

    /// **The tolerance the existing SSR caller states**, `1e-6`
    /// ([`em.rs`](../../../../src/ssr/cohort/em.rs)'s `EmCfg::tol`).
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
            crate::ng::parameter_estimation::ssr::StratumKey,
            crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
        >,
        calibration: &[ReadGroupCalibration],
        convergence_threshold: f64,
        scratch: &mut CallingScratch<StutterSubstitutionScratch>,
    ) -> LocusInference {
        let detail = tract_detail();
        let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
            .iter()
            .map(|observations| SsrSampleEvidence::new(observations, &detail))
            .collect();
        let repeat_counts = tract_repeat_counts();
        let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);
        let inbreeding = vec![outbred(); per_sample.len()];
        let parameters =
            uncontaminated_run(calibration, &inbreeding, strata, substitution, diploid());
        let config = CallingLoopConfig {
            convergence_threshold,
            ..CallingLoopConfig::DEFAULT
        }
        .validate()
        .expect("a threshold inside the configured range");
        shipped_arm().call_locus(&evidence, &parameters, tract_alleles(), &config, scratch)
    }

    /// **The same tract with the slippage round switched on** — `rounds` re-fit rounds at the
    /// given pull-backs, everything else the shipped configuration. What the re-fit fixtures
    /// call, and nothing else uses.
    fn call_tract_with_refit(
        observations_of_each_sample: &[&[SequenceObservation]],
        rounds: u32,
        direction_and_fall_off_pull_back: f64,
        level_pull_back: f64,
        scratch: &mut CallingScratch<StutterSubstitutionScratch>,
    ) -> LocusInference {
        let detail = tract_detail();
        let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
            .iter()
            .map(|observations| SsrSampleEvidence::new(observations, &detail))
            .collect();
        let repeat_counts = tract_repeat_counts();
        let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);
        let inbreeding = vec![outbred(); per_sample.len()];
        let strata = tract_strata();
        let substitution = tract_substitution_rates();
        let calibration = tract_libraries();
        let parameters =
            uncontaminated_run(&calibration, &inbreeding, &strata, &substitution, diploid());
        let config = CallingLoopConfig {
            slippage_refit: super::super::SlippageRefitConfig {
                max_rounds: rounds,
                direction_and_fall_off_pull_back_pseudocounts: direction_and_fall_off_pull_back,
                level_pull_back_slipped_reads: level_pull_back,
                round_convergence_threshold: super::super::DEFAULT_ROUND_CONVERGENCE_THRESHOLD,
            },
            ..CallingLoopConfig::DEFAULT
        }
        .validate()
        .expect("a built setting with legal pull-backs");
        shipped_arm().call_locus(&evidence, &parameters, tract_alleles(), &config, scratch)
    }

    /// The three numbers the tract was last scored under at `(read group 0, candidate 0)` —
    /// the one cell every re-fit fixture's reads land in.
    fn effective_cell_0_0(
        scratch: &CallingScratch<StutterSubstitutionScratch>,
    ) -> crate::ng::parameter_estimation::joint::ssr_fit::Slippage {
        scratch
            .tract_fits()
            .effective_slippage_of_cell(ReadGroupId(0), 0)
            .expect("the fixture's fit reaches this cell")
    }

    /// **A tract with planted heavy slippage pulls its level toward the reads' own rate and
    /// is held short of it by the 20-slipped-reads pull-back**, with the expected numbers
    /// computed by hand (spec §5.1).
    ///
    /// One sample: 30 reads at the reference's six repeats and 10 reads at five — one whole
    /// unit short, a length no candidate carries, so under the near-certain `0/0` call every
    /// one of them is a one-unit contraction slip. The frozen cell (read group 0, stratum 6)
    /// says `level 0.04, shorter_share 0.83, fall_off 0.25`.
    ///
    /// The hand arithmetic, with the genotype posterior ≈ 1 on `0/0` (thirty clean reads
    /// against ten one-step slips leave the alternatives tens of nats behind, so every weight
    /// below is within ~10⁻⁶ of the whole read count):
    ///
    /// - attributed weight ≈ 40, slipped ≈ 10, expected slips ≈ 40 × 0.04 = 1.6;
    /// - level multiplier = (10 + 20) / (1.6 + 20) = 1.3888…, so the level moves
    ///   0.04 → **0.0555…** — toward the reads' own rate of 10/40 = 0.25, and nowhere near
    ///   it: the pull-back holds it at just over a fifth of the way;
    /// - direction split = (10 + 50 × 0.83) / (10 + 50) = **0.8583…**;
    /// - every slip is one unit, so the raw fall-off estimate is 0 and the blend is
    ///   (0 + 50 × 0.25) / 60 = **0.2083…**.
    ///
    /// The second round attributes the same reads under the re-fitted table, lands on the
    /// same pooled counts — the expected-slips denominator reads the **frozen** level, so the
    /// multiplier does not compound — and stops for moving less than the threshold. Which is
    /// also the early stop of the rounds' own rule: one round adopted of the three allowed.
    #[test]
    fn planted_slippage_moves_the_level_toward_the_reads_but_the_pull_back_holds_it() {
        let sample = [
            tract_reads(TRACT_CANDIDATE_REPEATS[0], 30),
            tract_reads(TRACT_CANDIDATE_REPEATS[0] - 1, 10),
        ];
        let mut scratch = worker_scratch();
        let inference = call_tract_with_refit(&[&sample], 3, 50.0, 20.0, &mut scratch);

        assert_eq!(called_alleles(&inference, 0), vec![0, 0]);
        let refitted = effective_cell_0_0(&scratch);
        assert!(
            (refitted.level - 0.04 * (30.0 / 21.6)).abs() < 1e-3,
            "level {} against the hand-computed 0.05555…",
            refitted.level
        );
        assert!(
            refitted.level < 0.25,
            "the pull-back must hold the level short of the reads' own rate"
        );
        assert!((refitted.shorter_share - 51.5 / 60.0).abs() < 1e-3);
        assert!((refitted.fall_off - 12.5 / 60.0).abs() < 1e-3);

        // The rounds stopped early — one adopted of the three allowed — and the one adopted
        // round is visible as the one extra table build.
        assert_eq!(scratch.slippage_refit_rounds(), 1);
        assert_eq!(scratch.emission_cost().emission_builds, 2);
    }

    /// **Zero pull-back is the free setting: the numbers go to the reads' own rate**
    /// (HipSTR's behaviour, spec §5.1). Same reads as the planted fixture, so the hand
    /// arithmetic differs only in the pull-backs: the multiplier is 10 / 1.6 = 6.25, the
    /// level 0.04 × 6.25 = **0.25** — exactly 10 slipped reads in 40 — the direction split
    /// **1.0** (every slip was a contraction) and the fall-off **0** (every slip one unit).
    #[test]
    fn zero_pull_back_hands_the_numbers_to_the_reads() {
        let sample = [
            tract_reads(TRACT_CANDIDATE_REPEATS[0], 30),
            tract_reads(TRACT_CANDIDATE_REPEATS[0] - 1, 10),
        ];
        let mut scratch = worker_scratch();
        let inference = call_tract_with_refit(&[&sample], 3, 0.0, 0.0, &mut scratch);

        assert_eq!(called_alleles(&inference, 0), vec![0, 0]);
        let refitted = effective_cell_0_0(&scratch);
        assert!(
            (refitted.level - 0.25).abs() < 1e-3,
            "free level {} against the reads' own rate 0.25",
            refitted.level
        );
        assert!((refitted.shorter_share - 1.0).abs() < 1e-3);
        assert!(refitted.fall_off.abs() < 1e-3);
        assert_eq!(scratch.slippage_refit_rounds(), 1);
    }

    /// **A tract whose reads carry no slips stays at the stratum's values.** Four clean reads
    /// expect 4 × 0.04 = 0.16 slips, so the level multiplier is 20 / 20.16 = 0.9921 — and the
    /// shape has no slips to move at all. The largest move any cell's numbers make is the
    /// multiplier over the largest frozen level the gather holds (the second slippage group's
    /// seven-repeat cell, 0.105): 0.105 × 0.0079 ≈ 0.0008, **under the round threshold of
    /// 10⁻³** — so the very first round converges, its numbers are discarded (production's
    /// test-before-adopt), and the locus keeps the frozen values **exactly**: no round
    /// adopted, no table rebuilt.
    ///
    /// *Why four reads and not ten*: the stopping rule reads every cell, and at ten clean
    /// reads the multiplier is 20 / 20.4 = 0.9804, which moves that same 0.105 cell by
    /// 0.002 — over the threshold, so a round is adopted and the level eases off the frozen
    /// value. That is the behaviour production shares (its multiplier moves 0.0196 against
    /// the same 10⁻³ tolerance), not a defect; this fixture sits below the threshold so that
    /// *stays frozen* can be asserted bit for bit.
    #[test]
    fn a_tract_with_no_slips_stays_at_the_stratum_values() {
        let sample = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 4)];
        let mut scratch = worker_scratch();
        let inference = call_tract_with_refit(&[&sample], 3, 50.0, 20.0, &mut scratch);

        assert_eq!(called_alleles(&inference, 0), vec![0, 0]);
        let cell = effective_cell_0_0(&scratch);
        assert_eq!(cell.level, 0.04, "the frozen value, bit for bit");
        assert_eq!(cell.shorter_share, 0.83);
        assert_eq!(cell.fall_off, 0.25);
        assert_eq!(scratch.slippage_refit_rounds(), 0);
        assert_eq!(
            scratch.emission_cost().emission_builds,
            1,
            "a converging first round rebuilds nothing"
        );
    }

    /// **An ordinary site ignores the re-fit structurally** (spec §5.1's closing paragraph):
    /// the same SNP called with three rounds asked for and with none is the same call, the
    /// same single table build, and a round count of zero.
    #[test]
    fn an_ordinary_site_ignores_the_slippage_rounds() {
        let (alleles, _) = generic_locus(1);
        let carrier = [observation(1, 0, 8, 4, 4)];
        let per_sample = [called_sample(&carrier)];
        let evidence = LocusEvidence::generic(locus_region(), &per_sample);
        let libraries = [ReadGroupCalibration::defaulted()];
        let inbreeding = [outbred()];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let substitution = std::collections::BTreeMap::new();
        let parameters = FrozenParameters::uncontaminated(
            &libraries,
            &inbreeding,
            human_like_seed(),
            &strata,
            &substitution,
            diploid(),
        );
        let config = CallingLoopConfig {
            slippage_refit: super::super::SlippageRefitConfig {
                max_rounds: 3,
                ..super::super::SlippageRefitConfig::DEFAULT
            },
            ..CallingLoopConfig::DEFAULT
        }
        .validate()
        .expect("re-fit rounds are a built setting");
        let mut scratch = worker_scratch();
        let inference =
            shipped_arm().call_locus(&evidence, &parameters, alleles, &config, &mut scratch);

        assert_eq!(called_alleles(&inference, 0), vec![1, 1]);
        assert_eq!(scratch.slippage_refit_rounds(), 0);
        assert_eq!(scratch.emission_cost().emission_builds, 1);
    }

    /// The genotype a sample was called, as allele ids.
    fn called_alleles(inference: &LocusInference, sample: usize) -> Vec<u16> {
        inference.per_sample[sample]
            .genotype()
            .unwrap_or_else(|| panic!("sample {sample} was called missing"))
            .alleles()
            .iter()
            .map(|allele| allele.get())
            .collect()
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

    /// **The tract's prior seed is its stratum's fitted length spectrum, restricted to the
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
    fn the_tracts_prior_seed_is_its_fitted_length_spectrum_over_the_candidates() {
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
    /// is decided by is a length spectrum.**
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
        let calibration = [ReadGroupCalibration::defaulted(); TRACT_READ_GROUPS];
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
    /// The SNP/indel path has an oracle — the existing caller's own loop, handed the same
    /// likelihood table (`calling::loop_parity`). **The repeat-tract path does not**, and the
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
    ///   ([`em.rs`](../../../../src/ssr/cohort/em.rs)'s `run_pi_em`, whose `expected` **starts at**
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
    /// **five passes against three** on this tract — most of a run's extra work for the same
    /// three genotypes — and the two stopping points land **1.8 × 10⁻⁶ of a chromosome** apart,
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
        /// so the counts moved when that constant moved to 0.20 (the shipped rule took two
        /// passes at 0.05).
        const PASSES_AT_THE_TIGHT_TOLERANCE: u32 = 5;
        /// The shipped tolerance's count on the same tract.
        const PASSES_AT_THE_SHIPPED_TOLERANCE: u32 = 3;
        /// How far apart the two stopping points landed, per chromosome — see the doc above.
        const CHROMOSOMES_APART: f64 = 1.8e-6;

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
    /// A run whose repeat fit produced no length spectrum anywhere seeds every tract from the
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
            LengthSpectrumRung::StatedFlat,
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

    /// **A repeat tract's prior comes from its stratum's length spectrum, and the record says
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
        assert_eq!(rung_of(&fitted), LengthSpectrumRung::StratumsOwnFit);

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
            LengthSpectrumRung::StatedFlat,
            "the ladder always answers, and says from how far down"
        );
    }

    /// **The same, on a contaminated run** — because the contaminant seed adds three buffers
    /// that its uncontaminated twin — `a_second_tract_on_the_same_scratch_is_scored_on_its_own_parameters`,
    /// further down this module — leaves empty from beginning to end, so nothing there could
    /// show a clear that was dropped.
    ///
    /// **The defaulting tract runs first here too, and it does more work in this order.** The
    /// first tract's fit reaches no stratum at all, so its length spectrum is the flat one at a
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
    /// reached, in a run whose fit pooled that motif period's tracts.
    ///
    /// **It is the rung neither of the fixtures beside it can produce**, and each of them can
    /// be produced by a mistake: a lookup that ignored the pool answers `StatedFlat`, and one
    /// that read the pool where a stratum's own fit exists answers `PeriodsPooledTracts` at the
    /// fitted tract as well. So the three rungs are asserted at three fixtures rather than
    /// two (`population_diversity.md` §8, check 3).
    #[test]
    fn a_tract_whose_stratum_the_fit_missed_is_seeded_from_its_motif_periods_pool() {
        use crate::ng::parameter_estimation::joint::ssr_fit::PeriodLengthSpectrum;
        // A fit that reached no stratum of its own, but pooled every dinucleotide tract it saw.
        let pooled = StratumFits::over(&[], std::collections::BTreeMap::new())
            .with_period_length_spectra(std::collections::BTreeMap::from([(
                2,
                PeriodLengthSpectrum {
                    period: 2,
                    length_spectrum: vec![0.10, 0.20, 0.40, 0.20, 0.10],
                    concentration: 12.0,
                    tracts_fitted: 900,
                    strata_pooled: 7,
                    converged: true,
                },
            )]));

        let reads = [tract_reads(TRACT_CANDIDATE_REPEATS[0], 12)];
        let mut scratch = worker_scratch();
        let inference = call_tract(
            &[&reads],
            &pooled,
            &NO_SUBSTITUTION_RATES,
            &tract_libraries(),
            &mut scratch,
        );
        assert_eq!(rung_of(&inference), LengthSpectrumRung::PeriodsPooledTracts);
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
            crate::ng::parameter_estimation::ssr::StratumKey,
            crate::ng::parameter_estimation::Estimate<crate::ng::types::ErrorRate>,
        >,
        scratch: &mut CallingScratch<StutterSubstitutionScratch>,
    ) -> LocusInference {
        let detail = tract_detail();
        let per_sample: Vec<SsrSampleEvidence<'_>> = observations_of_each_sample
            .iter()
            .map(|observations| SsrSampleEvidence::new(observations, &detail))
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
        let evidence: [&[SequenceObservation]; 3] =
            [&reference_sample, &contaminated_sample, &carrier];

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

    /// **The contaminant seed's own numbers**, on a fitted length spectrum, over a candidate set
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
                bases: &clean,
                repeat_count: NonZeroU32::new(6).expect("six repeats"),
            },
            SsrCandidate {
                bases: &interrupted,
                repeat_count: NonZeroU32::new(4).expect("four whole repeats"),
            },
            SsrCandidate {
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
                length_spectrum: parameters
                    .ssr_length_spectrum_at(tract_detail().motif.ssr_period(), reference_repeats),
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
    /// term is the fit's length spectrum for this tract's stratum, frozen before calling
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

    /// A read group whose fit found the fraction below, measured on real evidence.
    fn contaminated_at(fraction: f64) -> ContaminationView {
        ContaminationView {
            fraction,
            markers_with_reads: 400,
            reads_on_markers: 1_000,
            source: ContaminationSource::ThisReadGroupsReads,
        }
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

    /// The frozen parameters for a contaminated run, over the default batching — one batch
    /// holding all of it, which is what a run that declared nothing gets.
    fn contaminated_run<'a>(
        calibration: &'a [ReadGroupCalibration],
        contamination: &'a [ContaminationView],
        batching: &'a SequencingBatches,
        inbreeding: &'a [InbreedingF],
        strata: &'a StratumFits,
    ) -> FrozenParameters<'a> {
        FrozenParameters::new(
            calibration,
            contamination,
            batching,
            inbreeding,
            human_like_seed(),
            strata,
            &NO_SUBSTITUTION_RATES,
            diploid(),
        )
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

        let reference_only =
            CandidateAlleles::new(Box::from(b"ACGT".as_slice()), LocusKind::Generic);
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
            (copies[0] - 0.007_744_253_786_433_548).abs() < 1e-12
                && (copies[1] - 1.992_255_746_213_566_4).abs() < 1e-12,
            "the one called sample's two chromosomes are nearly all alternative: {copies:?}"
        );
    }

    /// **A record may not claim its parameters were fitted when they were defaulted.**
    ///
    /// The fixture's calibration is `ReadGroupCalibration::defaulted()`, whose own provenance
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
            (outbred_quality - 32.323_853).abs() < 1e-3
                && (inbred_quality - 20.752_254).abs() < 1e-3,
            "the same reads under F = 0 and F = 0.9: {outbred_quality} and {inbred_quality}"
        );
        let copies = inference.cohort_expected_copies().copies();
        assert!(
            (copies[0] - 2.003_187_445_790_709_7).abs() < 1e-12,
            "scoring both rows at F = 0 gives 2.000584 instead: {copies:?}"
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
            witnessed_in_locus: WitnessedLocusPositions::from_half_open_runs([(0_u16, 1_u16)])
                .expect("one witnessed position"),
            read_group: ReadGroupId(0),
            bases: Box::from(b"A".as_slice()),
            num_reads: 50,
            q_sum: -100.0,
        }];
        let per_sample = [GenericLocusSample {
            evidence: GenericSampleEvidence::new(&rows, 0.0, &partials),
            genotype_must_be_missing: false,
        }];
        let evidence = LocusEvidence::generic(locus_region(), &per_sample);
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
        let calibration = [ReadGroupCalibration::defaulted()];
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
            (copies[0] - 1.019_039_125_169_859_6).abs() < 1e-12
                && (copies[1] - 0.980_960_874_830_140_3).abs() < 1e-12,
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
        let per_sample = [SsrSampleEvidence::new(&[], &detail)];
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
            SsrSampleEvidence::new(&[], &detail),
            SsrSampleEvidence::new(&[], &detail),
        ];
        let repeat_counts = tract_repeat_counts();
        let evidence = LocusEvidence::ssr(locus_region(), &per_sample, &detail, &repeat_counts);

        assert!(is_callable(&evidence, 0));
        assert!(is_callable(&evidence, 1));
    }

    // ────────────────────────────────────────────────────────────────────────────────
    // D2 — the two cost invariants: what a locus pays for, and what a pass does not
    // ────────────────────────────────────────────────────────────────────────────────

    /// One locus's evidence over `per_sample_observations` samples, where entry *s* says how
    /// many `(allele, read group)` rows sample *s* shows — **deliberately unequal**, because
    /// samples do not have equal observation counts and a fixture built so they do is the one
    /// shape that hides a three-way product (`spec/calling_em_loop.md` §13's test 5).
    fn rows_of_each_sample(per_sample_observations: &[usize]) -> Vec<Vec<GenericObservation>> {
        per_sample_observations
            .iter()
            .enumerate()
            .map(|(sample, observations)| {
                (0..*observations)
                    .map(|row| {
                        // Every row a different allele, and the read counts varied by sample so
                        // that no two rows are interchangeable. **What the fixture is for is the
                        // count, not the call** — and the calls it happens to produce run the
                        // other way from the obvious guess: the sample whose reads are spread
                        // over three alleles is called at 54.7 Phred and the one-row sample at
                        // 12.3, because spreading the reads is what makes a heterozygote's two
                        // alleles both visible.
                        observation(row as u16, 0, 4 + sample as u32, 2, 2)
                    })
                    .collect()
            })
            .collect()
    }

    /// **`candidates × Σ_s (observations in sample s)`, and `Σ_s` is the part that needs
    /// saying** — asserted as the formula rather than as a literal, at two loci that take
    /// different numbers of passes.
    ///
    /// Three samples showing one, two and three observations over three candidate alleles:
    /// the sum is `3 × (1 + 2 + 3)` = **18**. **The two wrong shapes this fixture separates,
    /// both measured on it:** charging the first row's count for every row gives
    /// `3 × 1 × 3` = **9**, and charging the locus's pooled total for every row gives
    /// `3 × 6 × 3` = **54**. A fixture whose three samples showed three observations each
    /// would report 27 under all three, which is exactly what the equal-count fixture beside
    /// this one cannot tell apart.
    ///
    /// **The same 18 at two passes and at four.** That is what "the table is paid for once"
    /// means: a build hoisted inside the frequency loop returns identical genotypes and costs
    /// the emission four times over, so nothing but the counter can tell the two apart.
    #[test]
    fn the_emission_count_is_candidates_times_the_sum_over_samples_at_every_pass_count() {
        let (alleles, _) = generic_locus(2);
        let rows = rows_of_each_sample(&[1, 2, 3]);
        let per_sample: Vec<GenericLocusSample<'_>> =
            rows.iter().map(|row| called_sample(row)).collect();
        let evidence = LocusEvidence::generic(locus_region(), &per_sample);
        let calibration = [ReadGroupCalibration::defaulted()];
        let inbreeding = vec![outbred(); 3];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let parameters = uncontaminated_run(
            &calibration,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
            diploid(),
        );

        let candidates = 3_u64;
        let observations: u64 = [1, 2, 3].iter().sum();
        let mut passes_seen = Vec::new();
        // **One scratch across both loci, which is the only shape a real worker has** — and
        // what makes the counter's own reset load-bearing. Measured with
        // `prepare_for_locus`'s reset deleted (2026-08-26): the second locus reports
        // `emission_builds: 2, emission_row_fills: 6, emission_evaluations: 36,
        // table_assemblies: 2, row_assemblies: 6` — every field doubled — and a fresh scratch
        // per locus hides it completely.
        let mut scratch = worker_scratch();
        for config in [capped_at(2), RunnableCallingLoopConfig::default()] {
            let inference = shipped_arm().call_locus(
                &evidence,
                &parameters,
                alleles.clone(),
                &config,
                &mut scratch,
            );
            passes_seen.push(inference.passes);
            assert_eq!(
                scratch.emission_cost(),
                EmissionCost {
                    emission_builds: 1,
                    emission_row_fills: 3,
                    emission_evaluations: candidates * observations,
                    table_assemblies: 1,
                    row_assemblies: 3,
                },
                "one emission build and one fold over three rows, {candidates} candidates × \
                 {observations} observations, at {} passes",
                inference.passes
            );
        }
        assert_eq!(
            passes_seen,
            vec![2, 4],
            "the two runs must differ in their pass count, or the invariant is untested"
        );
    }

    /// **No buffer of the worker's scratch moves or grows, however many passes the loop
    /// takes** — the loop's zero-allocation invariant, as far as a crate that forbids
    /// `unsafe` can observe it (`CallingScratch::buffer_fingerprints`).
    ///
    /// The same locus is called twice on **one** scratch, capped at two passes and then at
    /// the default, and every buffer comes back at the same address with the same length. A
    /// `Vec` that grew during a pass would have moved its bytes; a per-pass buffer added to
    /// the scratch would have changed a length.
    ///
    /// **What this cannot see is a temporary allocated and dropped inside a pass**, which
    /// leaves no trace anywhere — so that half is counted for real, in
    /// `tests/ng_calling_loop_allocation.rs`. It lives in a test binary of its own because a
    /// global allocator counts the whole process and the lib suite runs in parallel; this
    /// test is the cheap guard that runs on every build.
    #[test]
    fn no_buffer_of_the_scratch_moves_or_grows_however_many_passes_the_loop_takes() {
        let (alleles, _) = generic_locus(2);
        let rows = rows_of_each_sample(&[1, 2, 3]);
        let per_sample: Vec<GenericLocusSample<'_>> =
            rows.iter().map(|row| called_sample(row)).collect();
        let evidence = LocusEvidence::generic(locus_region(), &per_sample);
        let calibration = [ReadGroupCalibration::defaulted()];
        let inbreeding = vec![outbred(); 3];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let parameters = uncontaminated_run(
            &calibration,
            &inbreeding,
            &strata,
            &NO_SUBSTITUTION_RATES,
            diploid(),
        );
        let mut scratch = worker_scratch();

        let two_passes = shipped_arm().call_locus(
            &evidence,
            &parameters,
            alleles.clone(),
            &capped_at(2),
            &mut scratch,
        );
        let after_two = scratch.buffer_fingerprints();

        let four_passes = shipped_arm().call_locus(
            &evidence,
            &parameters,
            alleles,
            &RunnableCallingLoopConfig::default(),
            &mut scratch,
        );
        let after_four = scratch.buffer_fingerprints();

        assert_eq!(
            (two_passes.passes, four_passes.passes),
            (2, 4),
            "the two runs must differ in their pass count, or the invariant is untested"
        );
        assert_eq!(
            after_two, after_four,
            "every buffer is where it was and as long as it was: a Vec that grew during a \
             pass would have moved its bytes"
        );
    }

    /// **And the same on the contaminated path, which is the one that does per-pass work
    /// inside the loop.**
    ///
    /// The test above calls a run with no fraction fitted, so its loop reads a table nothing
    /// touches. With a fraction fitted the loop refills the batch copies, refills one sample's
    /// contaminant frequencies per row and assembles the whole table again, at every pass —
    /// three buffers written inside the loop and sized outside it, which is exactly what
    /// `buffer_fingerprints` claims to measure and what nothing measured until this fixture.
    #[test]
    fn no_buffer_moves_or_grows_on_a_contaminated_locus_either() {
        let (alleles, _) = generic_locus(2);
        let rows = rows_of_each_sample(&[1, 2, 3]);
        let per_sample: Vec<GenericLocusSample<'_>> =
            rows.iter().map(|row| called_sample(row)).collect();
        let evidence = LocusEvidence::generic(locus_region(), &per_sample);
        let calibration = [ReadGroupCalibration::defaulted()];
        let fractions = [contaminated_at(0.05)];
        let inbreeding = vec![outbred(); 3];
        let strata = StratumFits::over(&[], std::collections::BTreeMap::new());
        let batching = SequencingBatches::all_together_over(1, 3);
        let parameters =
            contaminated_run(&calibration, &fractions, &batching, &inbreeding, &strata);
        let mut scratch = worker_scratch();

        let fewer = shipped_arm().call_locus(
            &evidence,
            &parameters,
            alleles.clone(),
            &capped_at(2),
            &mut scratch,
        );
        let after_fewer = scratch.buffer_fingerprints();

        let more = shipped_arm().call_locus(
            &evidence,
            &parameters,
            alleles,
            &RunnableCallingLoopConfig::default(),
            &mut scratch,
        );
        let after_more = scratch.buffer_fingerprints();

        assert!(
            more.passes > fewer.passes,
            "the two runs must differ in their pass count, or the invariant is untested — \
             {} against {}",
            fewer.passes,
            more.passes
        );
        // **Sorted, because two of the buffers exchange pointers with every pass** — the M-step
        // swaps the cohort's expected copies with the previous pass's rather than copying
        // either, so the list's order depends on the parity of the pass count and its *contents*
        // are the claim. The uncontaminated fixture above compares unsorted only because both
        // its runs take an even number of passes.
        let mut after_fewer = after_fewer;
        let mut after_more = after_more;
        after_fewer.sort_unstable();
        after_more.sort_unstable();
        assert_eq!(
            after_fewer, after_more,
            "every buffer is where it was and as long as it was, over a loop that refills three \
             of them at every pass"
        );
    }

    /// **A wider locus does grow the buffers, and it is meant to** — the fingerprint is not
    /// a claim that the scratch never allocates, only that a *pass* does not.
    ///
    /// Without this, the test above passes against an implementation whose buffers never
    /// change because nothing ever prepares them.
    #[test]
    fn a_wider_locus_than_the_worker_has_met_does_grow_the_scratch() {
        let (narrow, _) = generic_locus(1);
        let (wide, _) = generic_locus(3);
        let rows = rows_of_each_sample(&[1]);
        let per_sample: Vec<GenericLocusSample<'_>> =
            rows.iter().map(|row| called_sample(row)).collect();
        let evidence = LocusEvidence::generic(locus_region(), &per_sample);
        let calibration = [ReadGroupCalibration::defaulted()];
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
            narrow,
            &RunnableCallingLoopConfig::default(),
            &mut scratch,
        );
        let after_narrow = scratch.buffer_fingerprints();

        let _ = shipped_arm().call_locus(
            &evidence,
            &parameters,
            wide,
            &RunnableCallingLoopConfig::default(),
            &mut scratch,
        );
        let after_wide = scratch.buffer_fingerprints();

        assert_ne!(
            after_narrow, after_wide,
            "a locus of four alleles is wider than one of two, so the per-genotype buffers \
             have to grow — a fingerprint that never moved would mean nothing was prepared"
        );
    }
}
