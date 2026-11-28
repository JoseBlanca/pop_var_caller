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
//! [`CallingScratch`] and is reused at every locus and
//! on every pass, which is why the arithmetic below takes a
//! [`SampleScoringBuffers`] rather than returning anything.

use crate::calling::GenotypeTableView;
use crate::calling::genotype_prior::seed_generic::{VariantClass, fill_locus_concentration};
use crate::calling::genotype_prior::{
    BundlePrior, CohortAlleleCopies, Concentration, GenotypePriorModel, PriorRow,
    SampleAlleleCopies, fill_sample_concentration, fill_ssr_bundle_seed, fill_ssr_seed,
};
use crate::calling::likelihood::generic::{
    assemble_genotype_log_likelihood_row, fill_generic_emissions,
};
use crate::calling::likelihood::ssr::genotype_log_likelihood_row;
use crate::calling::likelihood::ssr_emission::SsrEmissionModel;
use crate::calling::likelihood::{FrozenContamination, ReadGroupCalibrations};
use crate::calling::quality::{
    ArtifactTestCounts, SiteQualityPrior, reads_said_nothing, score_best_genotype,
    score_uncorrected_site_quality,
};
use crate::calling::{
    CallingScratch, CandidateAlleles, CohortSummingBuffers, ContaminationMixture,
    ExpectedAlleleCopies, FrozenParameters, GenericLocusSample, GenericSampleEvidence,
    GenotypeTable, LocusEvidence, LocusInference, ReadGroupParameters, RepeatTractProvenance,
    SampleGenotypeCall, SampleScoringBuffers, SsrLocusShape, SsrSampleEvidence,
    WrongReadChargeTable, fill_batch_allele_copies, fill_own_batch_contaminant_row,
    fill_wrong_read_charges,
};
use crate::float;
use crate::locus_generation::SsrBundleDetail;
use crate::parameter_estimation::Provenance;
use crate::parameter_estimation::repeat_strata::RepeatCount;
use crate::types::{AlleleId, Genotype, GenotypeRows, InbreedingF, LogProb, Ploidy};
use std::iter::repeat_n;

use super::repeat_tract_parameters::{
    TractPrior, TractScoringFits, bundle_candidates, tract_candidates,
};
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
/// (`d9e7b076:src/var_calling/posterior_engine.rs`, `EmStepPhase::FirstIteration`).
///
/// **Inherited from GATK and from production, and never measured here — soft.** The argument is
/// arithmetic about the prior's size, not a count of calls that moved. It should bite hardest
/// where the read likelihood is weakest against a 20-to-30 Phred prior, which is the tomato
/// panel's corner at three reads a position, and hardly at all at 300. Spec §12's question 7 is
/// the measurement.
#[derive(Clone, Copy)]
pub(crate) enum PassPrior<'a> {
    /// **No prior at all** — every candidate genotype equally likely, and no concentration
    /// built, so nothing reads the cohort's expected copies. The first pass at every locus.
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
///   sentinel makes this necessary: [`CallingScratch::prepare_for_locus`] fills that row with
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
// this module's, and the arm is the only thing outside it a run needs. Until the loop's driver
// wrote that arm the whole chain was dead in a non-test build and carried a `cfg_attr(not(test),
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
        *slot = float::exp(*slot - largest_score);
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
    //
    //    **A sample whose reads said nothing holds none**: its posterior is the prior alone, and
    //    counting it would feed the cohort's own frequency back in as if a sample had shown
    //    it — pulling the frequency toward the seed and moving every other sample's prior,
    //    with no read behind either. In the full model such a sample marginalises out exactly
    //    (`doc/devel/ng/spec/calling_em_loop.md` §7). Zero copies also keep its own
    //    leave-one-out prior right: the cohort less nothing is the cohort.
    sample_expected_copies.fill(0.0);
    if reads_said_nothing(genotype_likelihoods) {
        return;
    }
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
/// **Not handled here: it never reaches this table.** Spec §5.0 sets aside a sample whose own
/// reads earned an allele the cap then cut: it is scored against nothing and *contributes
/// nothing to this sum*, because its posterior would sit over the wrong set of alleles and
/// would pull the locus's frequencies toward the reference by exactly the samples carrying the
/// rarest alleles. The loop's driver ([`SummariseConditionLoop`]'s
/// [`call_locus`](LocusGenotyper::call_locus)) gives such a sample no scratch row at all, so
/// every row summed here belongs to a sample the loop scores. A row nobody wrote still holds
/// the `NaN` sentinel `prepare_for_locus` left, and the finiteness check below refuses it.
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
    // today, because `sample_count` is checked above — but it is exactly the shape a caller
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
    // **On the first pass at a locus only, and the difference matters to the loop's driver.**
    // The sentinel is written by `prepare_for_locus`, once per locus, and the per-sample rows
    // are deliberately *not* re-armed between passes, because `score_one_sample` reads each
    // sample's previous copies as its leave-one-out term. So from pass 2 a row this pass did
    // not write holds the **previous pass's finite, plausible value** and is summed as though
    // it were current — measured `[1.0, 3.0]` where pass 1 gave `[2.0, 2.0]` and one sample's
    // E-step did not run. Finite, non-negative, wrong, silent. The driver never skips a sample
    // mid-run — a set-aside sample has no row from the start — so no per-pass written mask is
    // needed; a loop that did skip one would need it.
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
/// identical trajectory (`d9e7b076:src/var_calling/posterior_engine.rs:2702`).
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
/// [`UNWRITTEN_SCRATCH_VALUE`](crate::calling::UNWRITTEN_SCRATCH_VALUE), which is `NaN`,
/// and [`advance_cohort_expected_copies`](crate::calling::CallingScratch::advance_cohort_expected_copies)
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
///   [`previous_cohort_expected_copies`](crate::calling::CallingScratch::previous_cohort_expected_copies)'s
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
/// advanced, such as the final pass.
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
/// travel into [`LocusInference`] unchanged.
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
/// must score them again. That is the final pass's job, and it is a pass rather than a
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
/// [`CallingScratch::inbreeding_coefficient_by_row`](crate::calling::CallingScratch)
/// refuses a scratch whose row map is not one entry per prepared row. It is the same check the
/// old slice argument carried, moved onto the scratch, where it fires at the **first read** of
/// the row map rather than in this function's preamble.
///
/// **The two directions fail differently, which is why the check covers both, and neither
/// fails the way this paragraph used to say** (measured with the check downgraded, in
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
    // counted as a pass — `passes` counts the passes that had a prior (§3).
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

        // **The samples whose reads said something**, since only they put copies in the sum:
        // measuring the movement against every row's chromosomes would loosen the threshold by
        // each sample with no reads here, and a cohort with one more such sample would stop at
        // a different frequency. A locus nobody covered has nothing moving, and settles on the
        // first pass at any count, so one stands in for none.
        let informative_rows = (0..row_count)
            .filter(|&row| !reads_said_nothing(scratch.sample_genotype_likelihoods(row)))
            .count()
            .max(1);
        if cohort_expected_copies_have_settled(
            scratch.previous_cohort_expected_copies(),
            scratch.cohort_expected_copies(),
            genotypes.ploidy(),
            informative_rows,
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
/// **The three read counts are summed as integers**, which is not a micro-optimisation but the
/// removal of a question: `f64` addition is not associative, so a cohort's read totals summed in
/// another order could differ in the last bits and the two binomial tests downstream would see
/// different inputs at a different worker count. Whole reads summed as `u64` cannot. **The other
/// five are fractional by construction** — a read composed across records carries a forward and
/// a placed-left *share* (`AlleleSupport::num_fwd`), and a genotype carrying one copy of two
/// expects half a sample's reads — so they are `f64`s summed in the run's fixed sample order,
/// the same order and the same reason as the M-step's (`doc/devel/ng/spec/calling_quality.md`
/// §9).
#[derive(Default)]
struct PooledArtifactCounts {
    reference_reads: u64,
    reference_forward_reads: f64,
    reference_placed_left_reads: f64,
    alternative_reads: u64,
    alternative_forward_reads: f64,
    alternative_placed_left_reads: f64,
    reference_low_quality_reads: u64,
    alternative_low_quality_reads: u64,
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
    /// (`d9e7b076:src/vcf/qual_refine.rs:92`).
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
                self.reference_forward_reads += observation.forward_reads;
                self.reference_placed_left_reads += observation.placed_left_reads;
                self.reference_low_quality_reads += u64::from(observation.low_quality_reads);
            } else if observation.allele == primary_alternative {
                self.alternative_reads += reads;
                self.alternative_forward_reads += observation.forward_reads;
                self.alternative_placed_left_reads += observation.placed_left_reads;
                self.alternative_low_quality_reads += u64::from(observation.low_quality_reads);
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
            reference_forward_reads: self.reference_forward_reads,
            reference_placed_left_reads: self.reference_placed_left_reads,
            alternative_reads: self.alternative_reads as f64,
            alternative_forward_reads: self.alternative_forward_reads,
            alternative_placed_left_reads: self.alternative_placed_left_reads,
            reference_low_quality_reads: self.reference_low_quality_reads as f64,
            alternative_low_quality_reads: self.alternative_low_quality_reads as f64,
            total_reads: self.total_reads as f64,
            genotype_expected_alternative_reads: self.genotype_expected_alternative_reads,
        }
    }
}

/// **Every allele's reads, pooled over the samples the locus is called on** — and whether any
/// alternative was reached at all, which is whether the artifact tests have anything to weigh.
///
/// `pooled_reads` is the scratch buffer this fills, one entry per allele of the locus; it
/// arrives holding the previous locus's totals and is zeroed here. Which alternative the tests
/// weigh is [`pick_primary_alternative`]'s, once the samples are called.
///
/// **`false` is two different situations and the caller need not tell them apart**: a locus
/// whose candidate table is the reference alone — 27.4% of built loci on the 63-accession
/// tomato panel and 27.3% on HG002 at 30×
/// ([`SelectionVerdict::Selected`](crate::calling::allele_candidates::SelectionVerdict))
/// — and a locus with alternatives that no read reached. Both leave the two tests with
/// nothing to weigh, and production returns its baseline unchanged in exactly these two
/// cases (`d9e7b076:src/vcf/qual_refine.rs:79`).
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
fn pool_allele_reads(per_sample: &[GenericLocusSample<'_>], pooled_reads: &mut [u64]) -> bool {
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
                observation.forward_reads <= f64::from(observation.num_reads)
                    && observation.placed_left_reads <= f64::from(observation.num_reads)
                    && observation.low_quality_reads <= observation.num_reads,
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
    // Allele 0 is the reference, which the tests weigh the alternative against.
    pooled_reads.iter().skip(1).any(|&reads| reads > 0)
}

/// **Which allele the artifact tests treat as *the* alternative**: of the alternatives some
/// sample was called carrying, the one the most reads reached — and only where no called
/// alternative was reached by a read, the most-read alternative of all.
///
/// **A called allele first, because the tests are about the record.** Both compare the reads
/// of one alternative against the reference, and the allele-balance test's expectation is the
/// copies of that alternative in the called genotypes. An alternative no sample carries is
/// not in the record's `ALT`, expects no reads, and penalises a site for evidence it does not
/// report: on 8 tomato accessions, SL4.0ch01:3493044 lost 126 of its 163 Phred to the strand
/// test on reads of an allele no sample was called with.
///
/// **Ties go to the alternative whose bases sort first**, never to the lower id. Ids follow
/// the order the merge first met the alleles in, which is the order the samples were given,
/// so an id tie-break made a site's quality depend on how the cohort was listed — 2 sites in
/// 3,200 on that slice, one of them written or not by it. Bases cannot tie.
///
/// `None` where no alternative was reached by a read; [`pool_allele_reads`] has already said so,
/// and this is not asked.
fn pick_primary_alternative(
    pooled_reads: &[u64],
    candidates: &CandidateAlleles,
    calls: &[SampleGenotypeCall],
) -> Option<AlleleId> {
    let carried_by_a_call = |allele: AlleleId| {
        calls.iter().any(|call| match call {
            SampleGenotypeCall::Called { genotype, .. } => genotype.alleles().contains(&allele),
            SampleGenotypeCall::Missing => false,
        })
    };
    let most_read = |eligible: &dyn Fn(AlleleId) -> bool| {
        let mut best: Option<(AlleleId, u64, &[u8])> = None;
        for (index, &reads) in pooled_reads.iter().enumerate().skip(1) {
            // PANIC-FREE: `pooled_reads` is one entry per allele of the genotype table, and
            // `GenotypeTable::build` refuses an `allele_count` above `MAX_ALLELE_COUNT` — so
            // the index is at most 65,535 and always fits a `u16`.
            let allele = AlleleId(u16::try_from(index).expect("an allele id fits a u16"));
            if reads == 0 || !eligible(allele) {
                continue;
            }
            let bases = candidates
                .bases_of(allele)
                .expect("the pooled reads are one entry per candidate allele");
            let better = best.is_none_or(|(_, best_reads, best_bases)| {
                reads > best_reads || (reads == best_reads && bases < best_bases)
            });
            if better {
                best = Some((allele, reads, bases));
            }
        }
        best.map(|(allele, _, _)| allele)
    };
    most_read(&carried_by_a_call).or_else(|| most_read(&|_| true))
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
/// and the cursor is the whole of the map back (spec §5.0, §9).
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
              locus carries, all of which the driver has in hand when it calls this; a bundle \
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
            reference_beside: _,
        } => Some(*per_sample),
        LocusEvidence::Ssr { .. } => None,
    };
    // **The rows the artifact summary pools from**, where some alternative was reached by a
    // read. Which alternative it weighs is chosen after the samples are called, because the
    // choice prefers an allele the calls carry.
    let artifact_rows = generic_samples
        .filter(|per_sample| pool_allele_reads(per_sample, scratch.pooled_allele_reads_mut()));

    // **The genotype table's ploidy, and there is no second source of it.** The table was
    // built for a `(ploidy, allele count)` shape and every genotype's copies come out of it;
    // taking `FrozenParameters`' would add a number nothing compares. It stays
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
        calls.push(SampleGenotypeCall::Called {
            genotype: mint_genotype(copies_of_each_allele),
            genotype_quality,
            // **Asked of the likelihood row, here, while this pass still has it.** The row is
            // per-sample scratch the next locus overwrites, so nothing downstream can recover
            // what the reads said as against what the prior assumed — which is the whole of
            // why the flag is minted at the point of scoring (`vcf::assemble`'s module note,
            // and `doc/devel/ng/spec/vcf_output.md` §7.1).
            reads_were_uninformative: reads_said_nothing(scratch.sample_genotype_likelihoods(row)),
            genotype_rows: GenotypeRows {
                ploidy,
                log_likelihoods: scratch
                    .sample_genotype_likelihoods(row)
                    .iter()
                    .map(|likelihood| likelihood.get())
                    .collect(),
                posteriors: scratch.posterior_row().into(),
            },
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

    // **The artifact summary, now that the calls exist to choose its alternative from.** Every
    // called sample's reads are pooled against that one allele — the samples the candidate
    // step set aside are `Missing` and contribute nothing, as `pool_allele_reads` left them out.
    //
    // **The rows and the allele are one value, because the second cannot exist without the
    // first.** Held as two `Option`s and re-paired at the use site, the combination that cannot
    // arise — an alternative allele with no rows behind it — would pool nothing and say
    // nothing, which is the silent-wrong-answer shape rather than the loud one.
    let artifact_pool = artifact_rows.and_then(|per_sample| {
        pick_primary_alternative(scratch.pooled_allele_reads_mut(), &candidates, &calls)
            .map(|primary| (per_sample, primary))
    });
    if let Some((locus_samples, primary)) = artifact_pool {
        for (locus_sample, call) in locus_samples.iter().zip(&calls) {
            if let SampleGenotypeCall::Called { genotype, .. } = call {
                let copies = genotype
                    .alleles()
                    .iter()
                    .filter(|&&allele| allele == primary);
                pooled.add_called_sample(
                    &locus_sample.evidence,
                    primary,
                    u32::try_from(copies.count()).expect("a genotype's copies fit a u32"),
                    ploidy,
                );
            }
        }
    }

    // **After the samples, and it could as well have been before them**: the fold reads the
    // genotype likelihood table, which the loop built and this pass only reads
    // (`doc/devel/ng/spec/calling_quality.md` §3.2).
    //
    // **The prior is the one this locus's genotypes were called under** (`calling_quality.md`
    // §8.1): the run's frequency spectrum at an ordinary locus, and at a repeat tract the seed
    // concentration its own allele length variation gave — still in the scratch, which nothing after
    // the per-locus fill writes. Until 2026-10-10 a tract read the run's spectrum here, about
    // 1 variant position in 1,000, and at one sample a heterozygote called at genotype quality
    // 20 came out of this fold at 0.6 and was not written.
    let prior = match evidence {
        LocusEvidence::Generic { .. } => SiteQualityPrior::from(parameters.prior_seed()),
        LocusEvidence::Ssr { .. } => {
            SiteQualityPrior::of_repeat_tract(scratch.seed_concentration())
        }
    };
    let site_quality =
        score_uncorrected_site_quality(scratch.site_quality_buffers_mut(), genotypes, prior);

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
/// **The loop's currency is the [`GenotypeIdx`](crate::calling::GenotypeIdx) and the
/// output's is this**, and the conversion happens once per called sample on the final pass
/// rather than anywhere in the loop: a `Genotype` owns a boxed slice, so minting one per
/// sample per *pass* would be an allocation on the hot path for a value every pass but the
/// last throws away (`doc/devel/ng/arch/calling_em_loop.md` §2).
fn mint_genotype(copies_of_each_allele: &[u32]) -> Genotype {
    let mut alleles = Vec::with_capacity(copies_of_each_allele.iter().sum::<u32>() as usize);
    for (allele, &copies) in copies_of_each_allele.iter().enumerate() {
        // PANIC-FREE: as at `pick_primary_alternative` — the row is one entry
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
            reference_beside: _,
        } => per_sample[run_sample].is_callable(),
        LocusEvidence::Ssr { .. } => true,
    }
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
/// provenance at all ([`SpectrumSeed`](crate::calling::genotype_prior::SpectrumSeed) is
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
            reference_beside: _,
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
            reference_beside: _,
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
/// *(One caller, the tract's table fill. The per-locus prior seed destructures the evidence
/// itself.)*
struct TractEvidence<'a> {
    /// One entry per sample of the run, in run order. **A tract rules no sample out**
    /// (`doc/devel/ng/spec/calling_em_loop.md` §5.0.1), so every sample has one.
    per_sample: &'a [SsrSampleEvidence<'a>],
    /// The repeat ground: a segment's motif and repeat counts, or a bundle's runs and parses.
    shape: SsrLocusShape<'a>,
}

/// **What the genotype prior believes about one repeat locus's lengths, looked up once per
/// locus**: a segment's [`TractPrior`] or a bundle's [`BundlePrior`].
///
/// Both the prior's own seed and the third term of the read-likelihood mixture read it, and they
/// must read the same one: the run reports which rung of the tract ladder answered, so a second
/// lookup keyed differently would make the run's own record false. [`str_prior_of`] is the one
/// lookup; [`SummariseConditionLoop::call_locus`] makes it before the seed and hands the same
/// value to the build of the tract's likelihood table.
enum StrLocusPrior<'a> {
    /// An STR segment: its stratum's allele length variation, keyed by the reference tract's
    /// repeat count.
    Segment(TractPrior<'a>),
    /// An STR bundle: each reference run's own count and its stratum's allele length variation.
    Bundle(BundlePrior<'a>),
}

/// **The repeat locus's prior, [`StrLocusPrior`]; `None` at an SNP/indel site**, which has none.
///
/// **The reference tract's repeat count is entry 0** of a segment's candidate counts, because the
/// reference allele is id 0 of every candidate table and these counts run parallel to it. It is
/// what the spectrum's offsets are measured from, and passing a candidate's instead would
/// re-centre the shape on that candidate. A bundle reads its reference runs the same way, from
/// entry 0 of its candidate parses ([`bundle_prior_of`]).
fn str_prior_of<'a>(
    evidence: &LocusEvidence<'_>,
    parameters: &'a FrozenParameters<'a>,
) -> Option<StrLocusPrior<'a>> {
    match evidence {
        LocusEvidence::Generic { .. } => None,
        LocusEvidence::Ssr {
            shape:
                SsrLocusShape::Segment {
                    detail,
                    candidate_repeat_counts,
                },
            ..
        } => {
            let reference_repeats = RepeatCount(candidate_repeat_counts[0].get());
            Some(StrLocusPrior::Segment(TractPrior {
                reference_repeats,
                allele_length_variation: parameters
                    .ssr_allele_length_variation_at(detail.motif.ssr_period(), reference_repeats),
            }))
        }
        LocusEvidence::Ssr {
            shape:
                SsrLocusShape::Bundle {
                    detail,
                    candidate_parses,
                },
            ..
        } => Some(StrLocusPrior::Bundle(bundle_prior_of(
            detail,
            candidate_parses,
            parameters,
        ))),
    }
}

/// **What an STR bundle's prior is read from, looked up once per locus**: each reference run's
/// own reference count and its stratum's allele length variation (`str_bundle_loci.md` §7.4). The
/// reference counts are the reference allele's run parse, entry 0 of the candidates.
fn bundle_prior_of<'a>(
    detail: &SsrBundleDetail,
    candidate_parses: &[crate::region_typing::run_parse::RunParse],
    parameters: &'a FrozenParameters<'a>,
) -> BundlePrior<'a> {
    let reference_run_counts: Vec<u32> = candidate_parses[0].copies().collect();
    let run_spectra = detail
        .runs
        .iter()
        .zip(&reference_run_counts)
        .map(|(run, &copies)| {
            parameters.ssr_allele_length_variation_at(run.motif().ssr_period(), RepeatCount(copies))
        })
        .collect();
    BundlePrior {
        reference_run_counts,
        run_spectra,
        stated_concentration: parameters.ssr_slippage_fits().stated_concentration(),
    }
}

/// The repeat-tract path's evidence, as [`generic_evidence_of`] is the SNP/indel path's — and
/// unreachable on that path for the same reason.
fn tract_evidence_of<'a>(evidence: &'a LocusEvidence<'a>) -> TractEvidence<'a> {
    match evidence {
        LocusEvidence::Ssr {
            region: _,
            per_sample,
            shape,
        } => TractEvidence {
            per_sample,
            shape: *shape,
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
/// true is [`EmissionCost`](crate::calling::EmissionCost), which this charges as it goes: a
/// version that recomputed the work on every pass would give identical genotypes, only slower,
/// so nothing but a counter can tell the two apart.
///
/// # What each path fills
///
/// **On the SNP/indel path, the emissions and nothing else.**
/// [`assemble_genotype_likelihood_table`] then folds them into a per-genotype row, once per
/// locus or once per pass depending on contamination. Three things are filled and none of them
/// reads a frequency: the locus's wrong-read charge table, which is a property of the candidate
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
/// which the loop rewrites at every pass; a tract's is the *fit's* allele length variation for this
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
///
/// `str_prior` is the locus's [`str_prior_of`], looked up once by the caller: `None` at an
/// SNP/indel site, and at a repeat locus the same value the prior's seed was built from.
fn fill_what_no_pass_recomputes<Model: SsrEmissionModel>(
    emission: &Model,
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    candidates: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    scratch: &mut CallingScratch<Model::Scratch>,
    str_prior: Option<&StrLocusPrior<'_>>,
) {
    match evidence {
        LocusEvidence::Generic { .. } => {
            fill_generic_locus_emissions(evidence, parameters, candidates, genotypes, scratch);
        }
        LocusEvidence::Ssr { .. } => {
            let prior =
                str_prior.expect("a repeat locus's caller looked up its prior with `str_prior_of`");
            fill_tract_likelihood_table(
                emission, evidence, parameters, candidates, genotypes, scratch, prior,
            );
        }
    }
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
    let reference_beside = match evidence {
        LocusEvidence::Generic {
            reference_beside, ..
        } => *reference_beside,
        LocusEvidence::Ssr { .. } => unreachable!("generic_evidence_of refused a repeat tract"),
    };

    // **Once per locus rather than once per sample**: how far an allele's own error mass is
    // spread across the locus's others depends on the candidate sequences and on nothing a
    // sample showed.
    fill_wrong_read_charges(candidates, genotypes, scratch.wrong_read_charges_mut());

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
            reference_beside,
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
/// per locus (`doc/devel/ng/spec/read_likelihoods.md` §4.4).
///
/// # What it allocates, and it is the repeat path's own cost
///
/// **Two vectors per tract**: the candidate table as the row sees it — one
/// [`SsrCandidate`](crate::calling::likelihood::ssr_emission::SsrCandidate) per candidate,
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
/// [`tract_candidates`].
///
/// On a contamination list that disagrees with the fit these parameters were gathered under —
/// by [`TractScoringFits::locus_parameters`](super::repeat_tract_parameters::TractScoringFits::locus_parameters).
///
/// On a `prior` of the other kind of repeat locus than `evidence` — a segment's prior with a
/// bundle's evidence, or the reverse. [`str_prior_of`] looks the prior up from the evidence, so
/// this is a caller bug.
#[allow(
    clippy::too_many_arguments,
    reason = "the same inputs as `fill_what_no_pass_recomputes`, which passes them on"
)]
fn fill_tract_likelihood_table<Model: SsrEmissionModel>(
    emission: &Model,
    evidence: &LocusEvidence<'_>,
    parameters: &FrozenParameters<'_>,
    candidates: &CandidateAlleles,
    genotypes: &GenotypeTableView<'_>,
    scratch: &mut CallingScratch<Model::Scratch>,
    prior: &StrLocusPrior<'_>,
) {
    let tract = tract_evidence_of(evidence);
    let scored_candidates = match (tract.shape, prior) {
        (
            SsrLocusShape::Segment {
                detail,
                candidate_repeat_counts,
            },
            StrLocusPrior::Segment(tract_prior),
        ) => {
            let scored = tract_candidates(candidates, candidate_repeat_counts);
            scratch.tract_fits_mut().gather_for_locus(
                &detail.motif,
                &scored,
                *tract_prior,
                parameters,
            );
            scored
        }
        (
            SsrLocusShape::Bundle {
                detail,
                candidate_parses,
            },
            StrLocusPrior::Bundle(bundle_prior),
        ) => {
            let scored = bundle_candidates(candidates, detail, candidate_parses);
            scratch
                .tract_fits_mut()
                .gather_for_bundle(detail, &scored, bundle_prior, parameters);
            scored
        }
        _ => unreachable!(
            "a repeat locus's prior is looked up from its own evidence by `str_prior_of`, so a \
             segment's evidence always arrives with a segment's prior and a bundle's with a \
             bundle's"
        ),
    };

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
    let mut contexts = buffers.fits.scoring_contexts(&scored_candidates);
    // **Each read group's own misread table after runs** (`str_context_error_fit.md` §6). The
    // contexts are read-group-major, the read group's index its identifier, as the gathering
    // filled them.
    for (cell, context) in contexts.iter_mut().enumerate() {
        let group = u32::try_from(cell / scored_candidates.len())
            .expect("a run has fewer read groups than a u32 can name");
        context.misread_after_runs = parameters.misread_table(crate::types::ReadGroupId(group));
    }
    let locus = buffers.fits.locus_parameters(
        &scored_candidates,
        &contexts,
        contamination_of_each_read_group,
    );
    let genotype_count = buffers.genotype_count;
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
    }
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
/// tract's third mixture term is the fit's allele length variation for its stratum, frozen before
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
            reference_beside: _,
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
    // **The frequencies every row shares, filled and checked once a pass**; each row then fills
    // and checks only its own batch's row, the one its subtraction changes. Lent out of the
    // scratch for the pass and given back after it.
    let shared = (frequencies == ContaminantFrequencies::TheLoopsOwnEstimate)
        .then(|| scratch.lend_shared_contaminant_frequencies().0);
    let checked = match (&frozen_contamination, &shared) {
        (Some(frozen), Some(shared)) => Some(frozen.check_frequencies(shared, allele_count)),
        _ => None,
    };

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
            // belongs to the run's report of what contamination it used, which does not report
            // it yet; it is bound here rather than dropped so that the next reader sees it is
            // owed rather than handled.
            let _own_batch_had_no_evidence = fill_own_batch_contaminant_row(
                buffers.batch_allele_copies,
                SampleAlleleCopies::new(buffers.own_expected_copies),
                own_batch,
                buffers.allele_count,
                buffers.own_contaminant_row,
            );
            *buffers.own_contaminant_batch = own_batch.get() as usize;
        }
        let buffers = scratch.generic_row_buffers_mut(row);
        let mixture = match (frozen_contamination, checked) {
            (Some(frozen), Some(checked)) => frozen.with_checked_frequencies(
                checked,
                parameters.batch_of_sample(run_sample),
                buffers.own_contaminant_row,
            ),
            _ => ContaminationMixture::uncontaminated(),
        };
        assemble_genotype_log_likelihood_row(
            &sample,
            genotypes,
            ReadGroupParameters::new(calibration, mixture),
            WrongReadChargeTable::over(buffers.wrong_read_charges, genotypes),
            buffers.row_scratch,
            buffers.genotype_likelihoods,
        );
    }
    if let Some(shared) = shared {
        scratch.return_shared_contaminant_frequencies(shared);
    }
}

/// **Summarise the cohort, then condition each sample on it** — arm A of the step-9 seam, and
/// the whole of `doc/devel/ng/spec/calling_em_loop.md` §2's pseudocode in one function.
///
/// One loop: the **frequency loop** ([`run_frequency_loop`]), over a likelihood table built
/// once per locus.
///
/// **Spec §2's pseudocode has two more loops around it, and neither is here.**
///
/// - **The discovery round (§4.1).** Building it showed that discovery's eligible set reads no
///   posterior, so discovery runs as a **pre-pass inside candidate selection** before this arm
///   is called (`DiscoveryMode::BeforeTheLoop`;
///   `doc/devel/ng/research/tract_genotype_accuracy_2026-09-03.md` §6.5), and the loop here
///   runs over a table already widened.
/// - **The slippage round (§5.1)**, which re-fitted a repeat tract's slippage numbers from its
///   own reads and rebuilt the table under them. It was built, measured and deleted (owner,
///   2026-11-28). Measured on the GIAB HG002 tandem-repeat benchmark at 30×, one sample, over
///   fitted per-stratum slippage: at the designed pull-backs (50 pseudo-counts on the direction
///   split and fall-off, 20 slipped reads on the level, at most 3 rounds) it changed the
///   genotype of **7 tracts in 6,993, net +1 correct** (4 spurious heterozygotes fixed, 2 real
///   ones collapsed, 1 call silenced). With no pull-back (HipSTR's setting) collapsed
///   heterozygotes rose from 28 to 72 and period-2+ accuracy fell 0.56 points
///   (`doc/devel/ng/research/tract_accuracy_program_report.md`, lever L3). It was never measured
///   at low depth or on a cohort; doing so means rebuilding it (the last commit that holds it
///   is `a9690a90`).
///
/// **The locus's emissions are computed exactly once**, before the frequency loop, whatever
/// the pass count — which is what makes the expensive half of the
/// likelihood's cost independent of it (§2, §8), and is what
/// [`EmissionCost`](crate::calling::EmissionCost) records. **The genotype-likelihood table
/// assembled from them is written once too where nothing is contaminated**, because with no
/// fraction fitted the assembled row reads no allele frequency; where something is, it is
/// written again at the head of every pass and once more against the frequencies the loop
/// settled on (`doc/devel/ng/spec/read_likelihoods.md` §3.6).
///
/// # A sample the candidate step ruled uncallable leaves before the first pass
///
/// It is given **no scratch row at all** (spec §5.0's ruling). The rows are the run's sample
/// order with the gaps closed up, so the M-step sums the rows there are, the convergence delta divides by the chromosomes those
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
    #[cfg(test)]
    fn name(&self) -> &'static str {
        "summarise the cohort, then condition each sample on it (arm A)"
    }

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
        // tract** it is not: the third term there is the *fit's* allele length variation for this
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
        // **allele length variation** — how that stratum's chromosomes are spread over tract lengths —
        // which is looked up per locus.
        //
        // **Both have a fallback ladder and only one of them reports which rung it used.** The
        // ordinary-site ladder is `population_diversity.md` §3.4's, as amended by
        // `ordinary_site_seed.md` §4, and its rung is resolved when the run's seed is projected
        // and is not carried onto a locus; the tract ladder's rungs are §4.4's three plus the
        // one `str_fit_to_calling.md` §3.2 adds between the top two, and the rung comes back on
        // the lookup and travels onto this locus's record.
        // **The repeat locus's prior, looked up once here** and read both by the seed below and
        // by every build of the tract's likelihood table ([`StrLocusPrior`] says why it must be
        // one lookup). `None` at an ordinary site.
        let str_prior = str_prior_of(evidence, parameters);
        let allele_length_variation_rung = match (evidence, &str_prior) {
            (LocusEvidence::Generic { .. }, _) => {
                let _ = fill_locus_concentration(
                    parameters.prior_seed(),
                    variant_class_of(&candidates),
                    candidates.len(),
                    scratch.seed_concentration_mut(),
                );
                None
            }
            (
                LocusEvidence::Ssr {
                    shape:
                        SsrLocusShape::Bundle {
                            candidate_parses, ..
                        },
                    ..
                },
                Some(StrLocusPrior::Bundle(prior)),
            ) => {
                // **The product of the runs' seeds** (`str_bundle_loci.md` §7.4): each run's
                // shape read off its own stratum at its own offset, held with the most variable
                // run's conviction — `fill_ssr_bundle_seed` says why.
                let run_counts: Vec<u32> = candidate_parses
                    .iter()
                    .flat_map(|parse| parse.copies())
                    .collect();
                let _ = fill_ssr_bundle_seed(&run_counts, prior, scratch.seed_concentration_mut());
                Some(prior.rung())
            }
            (
                LocusEvidence::Ssr {
                    shape:
                        SsrLocusShape::Segment {
                            candidate_repeat_counts,
                            ..
                        },
                    ..
                },
                Some(StrLocusPrior::Segment(prior)),
            ) => {
                // **The prior is keyed by the reference tract's own repeat count, entry 0 of
                // the candidates' counts** ([`str_prior_of`] says why that entry).
                //
                // **A run whose fit produced no allele length variation anywhere is called here, at
                // the ladder's bottom rung, rather than refused** — the owner's ruling of
                // 2026-08-27, taken where `population_diversity.md` §5's *refuse a tract in a
                // run carrying no tract parameters* met §4.4's *the ladder always answers*.
                // **Refusing would turn a whole class of runs into a hard failure for a
                // condition the output already states**: the rung on the locus's record says
                // `StatedFlat`, which is how a call resting on a stated constant stays
                // distinguishable from one resting on a measurement without re-running
                // anything. §5 and §6 of that spec carry the ruling.
                //
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
                    prior.reference_repeats.get(),
                    prior.allele_length_variation,
                    scratch.seed_concentration_mut(),
                );
                Some(prior.allele_length_variation.rung())
            }
            (LocusEvidence::Ssr { .. }, _) => unreachable!(
                "`str_prior_of` answers a repeat locus with the prior of its own kind of ground"
            ),
        };

        // **The expensive half, once per locus**: the emissions read no allele frequency, so
        // nothing the frequency loop does can move them.
        fill_what_no_pass_recomputes(
            &self.emission,
            evidence,
            parameters,
            &candidates,
            &genotypes,
            scratch,
            str_prior.as_ref(),
        );
        // **And the cheap half, once here and then once per pass where something is
        // contaminated.** This first one is what the prior-free initialisation pass reads, and
        // there is no estimate of the frequencies for it to score against — no pass has run —
        // so it scores the reads alone, which is that pass's whole job (§3).
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
        // **One more assembly against the frequencies the loop settled on**, so that the final
        // pass scores every sample against the same estimate the convergence test just accepted
        // rather than against the one the last pass started from. Where nothing is
        // contaminated the table reads no frequency and there is nothing to assemble again.
        //
        // **What it costs is worth naming**: the genotypes and the site quality then come from
        // a table one assembly *newer* than the one the convergence test looked at. At a locus
        // that settled, the two differ by less than the convergence threshold by construction;
        // at one that hit the pass cap, the difference is whatever the last pass moved — and
        // such a locus is emitted flagged, which is what that flag is for (§6).
        if let Some(reassembly) = reassembly.as_ref() {
            reassembly.assemble(
                &genotypes,
                ContaminantFrequencies::TheLoopsOwnEstimate,
                scratch,
            );
        }

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
        let repeat_tract = allele_length_variation_rung.map(|rung| {
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
mod tests;
