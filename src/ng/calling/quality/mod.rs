//! **How sure the caller is** — one number per sample, and one for the locus.
//!
//! Two of the three numbers `doc/devel/ng/spec/calling_quality.md` defines are computed
//! here, and they are here rather than downstream for the same reason: **their inputs stop
//! existing when the locus is released**. The third, the artifact correction that adjusts
//! the site quality for strand and allele-balance skew, is arithmetic on nine scalars and
//! belongs to the output stream (§3.4); this module builds the summary it will read and
//! nothing more.
//!
//! - **The genotype quality** — how much of a sample's genotype probability did *not* go to
//!   the genotype that won ([`score_best_genotype`], §4). The posterior row it reads is one
//!   reused buffer of the worker's scratch, so by the time the last sample has been scored
//!   the earlier samples' posteriors are gone. It is taken during the loop's final pass,
//!   sample by sample, as they are scored.
//! - **The site quality, before its artifact correction** — *given every sample's reads, how
//!   unlikely is it that the cohort carries no copy of any non-reference allele at all?*
//!   ([`score_uncorrected_site_quality`], §5). It reads the whole `samples × genotypes` likelihood
//!   table, which is per-worker scratch overwritten at the next locus, and its fold is
//!   quadratic in cohort size — so computing it downstream would mean both carrying half a
//!   megabyte per locus in flight at 3,000 samples and putting a quadratic computation on
//!   the run's one serial thread (§3.2).
//!
//! **The site quality is not "how many samples look variant", and the difference has a
//! failure attached.** Multiplying each sample's probability of being homozygous reference
//! is not a normalised posterior: every additional reference-looking sample multiplies in
//! another factor below one, so that quality *grows with cohort size at a locus nobody
//! carries*. The marginal this module computes does not — adding a reference-looking sample
//! to a sparse-variant prior adds almost no evidence either way (§5.1).
//!
//! **Nothing here allocates.** The four buffers the fold needs belong to the worker's
//! [`CallingScratch`](crate::ng::calling::CallingScratch) beside the loop's own, and are
//! sized once per locus.

pub mod artifact_correction;

use crate::genetics::{MIN_ALT_CONCENTRATION, lgamma};
use crate::ng::calling::genotype_prior::SpectrumSeed;
use crate::ng::calling::genotype_prior::dirichlet_multinomial::log_sum_exp_2;
use crate::ng::calling::{GenotypeIdx, GenotypeTableView};
use crate::ng::types::{AlleleId, DomainError, LogProb, Phred};

/// The most a per-sample genotype quality is allowed to reach.
///
/// **A convention, not a measurement.** GATK and bcftools both cap here, so a downstream
/// tool reading ng's output sees the range it expects. Nothing in this repository has
/// measured what a cap of 99 costs or buys, and marking it soft is the point of saying so
/// (`doc/devel/ng/spec/calling_quality.md` §4).
pub const MAX_GENOTYPE_QUALITY: f32 = 99.0;

/// The most a site quality is allowed to reach — and the answer to [`Phred`]'s refusal of
/// an infinite value.
///
/// [`Phred::from_log_prob`] returns [`DomainError::PhredInfinite`](crate::ng::types::DomainError::PhredInfinite)
/// for a log-probability of `−∞`, deliberately, and its own documentation says the
/// consumer's answer is *"cap at its own ceiling and carry on"*. **This module is that
/// consumer and this is the ceiling** — production's `QUAL_MAX`, inherited and soft
/// (`doc/devel/ng/spec/calling_quality.md` §5.3).
///
/// **A `NaN` is not capped.** An infinite quality is a real answer at a locus whose reads
/// exclude the reference outright; a `NaN` is a bug in the arithmetic above and must
/// surface as one, which [`Phred::try_new`] already does.
pub const MAX_SITE_QUALITY: f32 = 9999.0;

/// The pooled read counts the two artifact tests read, and the one number the called
/// genotypes contribute — **nine scalars, whatever the cohort size**.
///
/// The tests themselves are the output stream's (`doc/devel/ng/spec/calling_quality.md`
/// §3.4). What has to happen in the worker is the *summing*: every quantity below is pooled
/// across the samples the locus was called on, and the per-sample evidence it is pooled from
/// is released with the locus. Nine numbers cross the boundary instead of a cohort-shaped
/// table, which is what makes the downstream stage affordable at three thousand samples
/// (§3.3).
///
/// **Counts are `f64` rather than integers**, because the two binomial tests consume them as
/// rates and production's do the same; carrying them as integers here would move the
/// conversion to the reader and gain nothing.
#[derive(Clone, Copy, PartialEq, Debug)]
pub struct ArtifactTestCounts {
    /// Which allele the tests treat as *the* alternative: the non-reference allele the most
    /// reads across the cohort reached.
    ///
    /// **One alternative, not all of them**, because both tests are two-sample comparisons
    /// against the reference and a triallelic locus has one quality, not three (§5.1).
    pub primary_alternative: AlleleId,
    /// How many reads across the cohort showed the reference allele.
    pub reference_reads: f64,
    /// Of those, how many were on the forward strand.
    pub reference_forward_reads: f64,
    /// Of those, how many started strictly left of the record they were seen at —
    /// freebayes' `placedLeft`, the read-position term.
    pub reference_placed_left_reads: f64,
    /// How many reads across the cohort showed [`Self::primary_alternative`].
    pub alternative_reads: f64,
    /// Of those, how many were on the forward strand.
    pub alternative_forward_reads: f64,
    /// Of those, how many were placed left.
    pub alternative_placed_left_reads: f64,
    /// Every allele's reads, summed over the samples the locus was called on — the
    /// denominator the allele-balance test compares against.
    pub total_reads: f64,
    /// How many alternative-allele reads the **called genotypes** lead you to expect:
    /// `Σ_s (copies of the primary alternative in sample s's call ÷ ploidy) × sample s's
    /// depth`.
    ///
    /// **The one entry the reads alone cannot give**, which is why the summary is built in
    /// the final pass rather than at the input edge: it needs the calls, and the calls are
    /// what the final pass produces.
    pub genotype_expected_alternative_reads: f64,
}

/// **How sure the caller is of one sample's genotype** — the winning genotype, and how much
/// of that sample's probability did not go to it.
///
/// ```text
/// GQ = min( 99, −10·log₁₀( 1 − p_best ) )
/// ```
///
/// Returns the winner's index alongside the quality because the two come from one walk of
/// the row and the caller needs both: the index mints the owned genotype, the quality goes
/// beside it (`doc/devel/ng/spec/calling_quality.md` §4).
///
/// # Three details, all inherited, and each is a trap without its reason
///
/// **`p_best` is nudged below one before the logarithm.** A sample whose reads make every
/// genotype but one impossible produces a posterior of exactly `1.0`, and `log₁₀(0)` is
/// `−∞`. Production clamps to one unit in the last place and so does this; the clamp is why
/// the result is always a real quality and never a [`DomainError`](crate::ng::types::DomainError).
///
/// **Ties go to the lower genotype index**, because the fold keeps the first *strict*
/// maximum. The genotype table's order is fixed, so this is deterministic — and it has to
/// be: two equally probable genotypes are exactly the case where an implementation that
/// kept the last maximum would give a different answer on a different machine.
///
/// **The cap is a convention** ([`MAX_GENOTYPE_QUALITY`]), not a measurement.
///
/// # Panics
///
/// Held in release, because a caller bug in this module is an assertion rather than a
/// `Result` (`doc/devel/ng/spec/calling_em_loop.md` §8):
///
/// - **the row must name at least one genotype.** A fold over an empty row would return the
///   index it started from and a quality computed from `−∞`, which is a call at a locus with
///   nothing to call.
/// - **the row must be a distribution: every entry finite, and the total one.** The E-step
///   normalises it, so this is a statement about that normalisation rather than about the
///   caller — and **the total is what the check is on, because the winner alone cannot see a
///   `NaN`.** No comparison against a `NaN` is true, so a `NaN` never wins the fold; a row of
///   `[0.7, NaN, 0.1]` returns genotype 0 at 5.229 Phred, a perfectly ordinary-looking call
///   made against a row one of whose genotypes has no probability at all. Summing catches it,
///   because `NaN` survives addition — and it catches a row that does not sum to one, which
///   is the other way a posterior can arrive wrong.
#[must_use]
pub(crate) fn score_best_genotype(posterior_row: &[f64]) -> (GenotypeIdx, Phred) {
    assert!(
        !posterior_row.is_empty(),
        "a locus has at least one candidate genotype, so a posterior row of none is a \
         scratch that was never prepared for this locus"
    );
    // The first *strict* maximum, so ties go to the lower index — and the total alongside it,
    // in the same walk, because that is the only thing that sees a `NaN` in a genotype that
    // did not win. `fold` rather than `max_by`, because `f64` has no total order and
    // `partial_cmp` on a `NaN` would have to be unwrapped.
    let (winner, best, total) = posterior_row.iter().enumerate().fold(
        (0_usize, f64::NEG_INFINITY, 0.0),
        |(best_index, best_probability, total), (index, &probability)| {
            if probability > best_probability {
                (index, probability, total + probability)
            } else {
                (best_index, best_probability, total + probability)
            }
        },
    );
    assert!(
        (total - 1.0).abs() < 1e-9,
        "a sample's posterior over genotypes is a distribution, so its {} entries must total \
         one: they total {total}. A NaN in any genotype reaches here and nowhere else — it \
         never wins the comparison that picks the best, since no comparison against a NaN is \
         true",
        posterior_row.len()
    );
    assert!(
        best.is_finite() && (0.0..=1.0).contains(&best),
        "the winning genotype's posterior probability is a probability, so it must be finite \
         and within [0, 1]: got {best} over {} genotypes",
        posterior_row.len()
    );

    // One unit in the last place below one, so the logarithm below is finite even where the
    // reads made every other genotype impossible.
    let best_below_one = best.min(1.0 - f64::EPSILON);
    let quality = -10.0 * (1.0 - best_below_one).log10();
    let capped = (quality as f32).clamp(0.0, MAX_GENOTYPE_QUALITY);
    (
        GenotypeIdx(u32::try_from(winner).expect("a genotype index fits a u32")),
        Phred::try_new(capped).expect("a value clamped into [0, MAX_GENOTYPE_QUALITY]"),
    )
}

/// The four buffers the site quality's fold reads and writes, borrowed from one
/// [`CallingScratch`](crate::ng::calling::CallingScratch) in one call.
///
/// **A bundle for the same reason the loop's two halves have one**: every field below is a
/// different field of one scratch, and reaching for them an accessor at a time does not
/// compile. Nothing here can alias.
/// **`pub` for one reader outside the crate: `benches/ng_site_quality_perf.rs`**, which times
/// the fold this bundle feeds across the cohort sizes the caller commits to
/// (`doc/devel/ng/spec/calling_quality.md` §13's Q3). The fields stay `pub(crate)` — a
/// benchmark passes the bundle on, it does not build one.
pub struct SiteQualityBuffers<'a> {
    /// How many samples the locus was called on.
    pub(crate) sample_count: usize,
    /// `samples × genotypes`, sample-major — the table the loop already built and never
    /// rebuilt. Read once, in the collapse.
    pub(crate) genotype_likelihoods: &'a [LogProb],
    /// `samples × (ploidy + 1)`, sample-major: each sample's log-likelihood of carrying
    /// exactly `c` non-reference copies, for `c` in `0..=ploidy`.
    pub(crate) copy_count_log_likelihoods: &'a mut [f64],
    /// The count axis the fold alternates between, `ploidy + samples × ploidy + 1` long.
    /// **Padded by `ploidy` on the left** so a copy-count tap can read `[padding − c]` without a
    /// bounds test in the inner loop.
    pub(crate) allele_count_distribution: &'a mut [f64],
    /// The other half of the alternation, the same length.
    pub(crate) allele_count_distribution_next: &'a mut [f64],
    /// `samples × ploidy + 1`: the fold's result back in the log domain, then the
    /// unnormalised log-posterior once the prior has been applied.
    pub(crate) log_allele_count_distribution: &'a mut [f64],
}

/// **How unlikely is it that the cohort carries no copy of any non-reference allele here?**
/// — the site quality, before the artifact correction that the output stream applies.
///
/// ```text
/// QUAL = −10·log₁₀ P( cohort allele count = 0 | every sample's reads )
/// ```
///
/// # What it is a statement about, and what it is not
///
/// **Not "is this sample variant" summed up.** The obvious alternative — multiply every
/// sample's probability of being homozygous reference — is not a normalised posterior, and
/// each reference-looking sample multiplies in another factor below one, so *that* quality
/// climbs with cohort size at a locus nobody carries. This one is the marginal posterior on
/// the cohort's allele count, so adding a reference-looking sample to a sparse-variant prior
/// adds almost no evidence either way and the number stays bounded by what the few
/// non-reference samples justify. Production replaced the first formula with this one and its
/// comment records the same reason; it is also what GATK means by `P(AC = 0 | data)`
/// (`doc/devel/ng/spec/calling_quality.md` §5.1).
///
/// **One quality over the union of the alternatives.** A triallelic locus does not get three
/// — *is there a variant here* is what a VCF's `QUAL` column asks, where `AC` and `AF` are
/// the per-allele questions and stay per-allele.
///
/// # Four steps (§5.2)
///
/// 1. **Collapse** each sample's row over genotypes into a row over *how many non-reference
///    copies that genotype carries* — `ploidy + 1` entries, by log-sum-exp over the genotypes
///    that share a count. The only step that reads the whole table.
/// 2. **Fold** the samples one at a time into a running distribution over the cohort's total.
/// 3. **Apply the prior** on that total (§5.4).
/// 4. **Normalise** over every possible total and read off the entry at zero.
///
/// # The prior is the run's own fitted spectrum, and that is where ng differs
///
/// The collapse to (reference, any-non-reference) turns the per-allele Dirichlet into a Beta
/// on the non-reference frequency, which induces a Beta-Binomial on the count. It needs two
/// concentrations. **Production's are two constants copied from GATK, each carrying "revisit
/// against the cohort calibration set" in its doc comment, and nobody did.** ng already holds
/// the same two numbers fitted from the run's own cohort — [`SpectrumSeed`] — so it uses
/// those (§5.4).
///
/// **It is not a free swap.** The prior sets how much read evidence a locus must produce
/// before its quality climbs off zero, and the two disagree by more as the cohort grows: at
/// 3,000 diploid samples the fitted seed asks for about 8 Phred more than production's
/// constants. That is the property worth having — on a panel ten times as polymorphic a locus
/// being variant is genuinely more likely before any read is looked at, and two fixed
/// constants cannot know that.
///
/// **The alternative concentration is floored at
/// [`MIN_ALT_CONCENTRATION`](crate::genetics::MIN_ALT_CONCENTRATION).** [`SpectrumSeed`]
/// admits exactly zero — a fully invariant cohort is a real answer — and `ln Γ(0)` is `+∞`,
/// which would make every term of the prior `NaN`. The floor sits far below any real
/// diversity, so it never moves a fitted estimate.
///
/// # Two numerical devices, and dropping either one silently ruins the answer
///
/// **The fold runs in the linear domain.** The mathematically identical log-domain version
/// spends its whole quadratic inner loop in `exp`/`ln`; production measured that at **88% of
/// its own path's time at 200 samples** before rewriting it as a plain multiply-add. Each
/// sample's copy-count row is divided by its own maximum, the folded result by *its* maximum, and
/// both divisions accumulate into a running log scale — which keeps every live value in
/// `(0, 1]` and leaves one `exp` per copy count and one `ln` per sample.
///
/// **The entry at zero is tracked separately, in logs.** It is the one the whole calculation
/// exists to read and the first to underflow: `log P(count = 0)` is `Σ_s (that sample's
/// log-likelihood of zero non-reference copies)`, a running sum, and at a confident cohort it
/// goes far below anything the rescaled linear buffer can hold. Read it back off the linear
/// array instead and **a strongly supported variant reports the ceiling** — which is why the
/// exact log value overrides the array's entry.
///
/// **Measured on this module's own fixture, samples each favouring the heterozygote by 20
/// nats:** at 20 samples the override changes nothing (1694.18 either way), at 50 it is the
/// difference between 4295.97 and the 9999 ceiling, and at 80 between 6899.70 and 9999. So it
/// starts to matter between 20 and 50 samples — well inside the range this caller commits to,
/// and below the 63 of the tomato panel.
///
/// # Panics
///
/// Held in release (`doc/devel/ng/spec/calling_em_loop.md` §8). Every one is a shape the
/// caller got wrong, and every one of them would otherwise produce a quality rather than a
/// failure: the sample count must be at least one; the likelihood table must be
/// `samples × genotypes`; each of the four buffers must be the length its field documents;
/// and the quality that comes out must not be a `NaN`, which is the only outcome of this
/// arithmetic that is a bug in it rather than an answer.
#[must_use]
pub fn score_uncorrected_site_quality(
    buffers: SiteQualityBuffers<'_>,
    genotypes: &GenotypeTableView<'_>,
    seed: SpectrumSeed,
) -> Phred {
    let SiteQualityBuffers {
        sample_count,
        genotype_likelihoods,
        copy_count_log_likelihoods,
        allele_count_distribution,
        allele_count_distribution_next,
        log_allele_count_distribution,
    } = buffers;

    let ploidy = usize::from(genotypes.ploidy().get());
    let genotype_count = genotypes.genotype_count();
    let allele_count = genotypes.allele_count();
    let copy_counts_per_sample = ploidy + 1;

    assert!(
        sample_count > 0,
        "a cohort has at least one sample, so a site quality over none is a run whose \
         sample order went missing"
    );
    let largest_count = sample_count
        .checked_mul(ploidy)
        .expect("a cohort of this many chromosomes cannot be indexed");
    assert_eq!(
        genotype_likelihoods.len(),
        sample_count * genotype_count,
        "the genotype likelihoods are samples × genotypes, sample-major: {sample_count} \
         samples over {genotype_count} genotypes is {} entries and the table holds {}",
        sample_count * genotype_count,
        genotype_likelihoods.len()
    );
    assert_eq!(
        copy_count_log_likelihoods.len(),
        sample_count * copy_counts_per_sample,
        "the copy-count log-likelihoods are samples × (ploidy + 1): {sample_count} samples at ploidy \
         {ploidy} is {} entries and the buffer holds {}",
        sample_count * copy_counts_per_sample,
        copy_count_log_likelihoods.len()
    );
    assert_eq!(
        log_allele_count_distribution.len(),
        largest_count + 1,
        "the count axis runs from 0 to ploidy × samples inclusive, which is {} entries; the \
         log buffer holds {}",
        largest_count + 1,
        log_allele_count_distribution.len()
    );
    for (which, axis) in [
        ("current", &*allele_count_distribution),
        ("next", &*allele_count_distribution_next),
    ] {
        assert_eq!(
            axis.len(),
            ploidy + largest_count + 1,
            "the {which} count axis is padded by ploidy on the left so a copy-count tap needs no \
             bounds test, which is {} entries; it holds {}",
            ploidy + largest_count + 1,
            axis.len()
        );
    }

    // 1. Collapse: genotypes → how many non-reference copies they carry.
    copy_count_log_likelihoods.fill(f64::NEG_INFINITY);
    let copies_per_genotype = genotypes.genotype_allele_counts();
    for (sample, copy_counts) in copy_count_log_likelihoods
        .chunks_exact_mut(copy_counts_per_sample)
        .enumerate()
    {
        let row = &genotype_likelihoods[sample * genotype_count..(sample + 1) * genotype_count];
        for (genotype, likelihood) in row.iter().enumerate() {
            // Allele 0 is the reference, always, so the non-reference copies are the rest.
            let reference_copies = copies_per_genotype[genotype * allele_count] as usize;
            let non_reference_copies = ploidy - reference_copies;
            let slot = &mut copy_counts[non_reference_copies];
            *slot = log_sum_exp_2(*slot, likelihood.get());
        }
        // **A sample with no possible genotype is caught here, where its cause has a name.**
        // Every entry of its collapsed row is then `−∞`, and the fold turns that into a
        // `−∞ − −∞` — a `NaN` that reaches the end of this function and fails
        // `Phred::from_log_prob` with a message about the *normalisation*, which is two steps
        // downstream of the mistake and sends the reader to the wrong place.
        assert!(
            copy_counts
                .iter()
                .any(|likelihood| *likelihood > f64::NEG_INFINITY),
            "sample {sample}'s reads make every one of the {genotype_count} candidate \
             genotypes impossible, so this locus has nothing to say about it; a likelihood \
             row that is entirely -inf reached the loop rather than being set aside"
        );
    }

    // 2. Fold the samples into the cohort's count distribution.
    fold_samples_into_allele_counts(
        ploidy,
        copy_count_log_likelihoods,
        allele_count_distribution,
        allele_count_distribution_next,
        log_allele_count_distribution,
    );

    // 3. The Beta-Binomial prior on the cohort's count, from the run's fitted spectrum.
    let alpha_reference = seed.alpha_ref();
    let alpha_alternative = seed.alpha_alt_total().max(MIN_ALT_CONCENTRATION);
    let chromosomes = largest_count as f64;
    let log_beta_normaliser = lgamma(alpha_alternative) + lgamma(alpha_reference)
        - lgamma(alpha_alternative + alpha_reference);
    let log_chromosomes_factorial = lgamma(chromosomes + 1.0);
    let log_denominator = lgamma(alpha_alternative + alpha_reference + chromosomes);
    for (count, slot) in log_allele_count_distribution.iter_mut().enumerate() {
        let count = count as f64;
        // How many ways a total of `count` can be dealt across the chromosomes, and the
        // Beta term for that split — the two halves of the Beta-Binomial.
        let log_ways =
            log_chromosomes_factorial - lgamma(count + 1.0) - lgamma(chromosomes - count + 1.0);
        let log_split =
            lgamma(alpha_alternative + count) + lgamma(alpha_reference + chromosomes - count);
        *slot += log_ways + log_split - log_denominator - log_beta_normaliser;
    }

    // 4. Normalise, and read the entry the whole calculation exists for.
    let log_total = log_sum_exp_over(log_allele_count_distribution);
    let log_probability_of_no_variant = log_allele_count_distribution[0] - log_total;
    match Phred::from_log_prob(LogProb(log_probability_of_no_variant)) {
        Ok(quality) if quality.get() <= MAX_SITE_QUALITY => quality,
        // A locus whose reads exclude the reference outright gives a probability of zero
        // here, which the Phred scale has no number for. Capping is this module's answer,
        // named once in `MAX_SITE_QUALITY` (§5.3).
        Ok(_) | Err(DomainError::PhredInfinite) => {
            Phred::try_new(MAX_SITE_QUALITY).expect("a finite non-negative ceiling")
        }
        Err(error) => panic!(
            "the site quality came out of a probability the Phred scale refuses: {error}. \
             The log-probability of no variant was {log_probability_of_no_variant} over {} \
             counts — above zero means the normalisation did not normalise, and a NaN means \
             the fold produced one",
            log_allele_count_distribution.len()
        ),
    }
}

/// Fold every sample's copy-count row into the cohort's allele-count distribution, writing the result
/// into `log_allele_count_distribution` in the log domain.
///
/// **Entry zero is not read back off the linear buffer.** It is accumulated exactly, in logs,
/// as `Σ_s (that sample's log-likelihood of zero non-reference copies)`, and overwrites
/// whatever the rescaled fold left there — see [`score_uncorrected_site_quality`] for why. See [`score_uncorrected_site_quality`]'s own comment for
/// why the fold is linear and why entry zero is tracked apart from it.
// **The two axis buffers share one lifetime**, which is what lets the loop below swap them.
// Two independent borrows cannot be swapped: `mem::swap` needs both sides to be the same
// type, and `&'a mut [f64]` and `&'b mut [f64]` are not until `'a` and `'b` are one.
fn fold_samples_into_allele_counts<'a>(
    ploidy: usize,
    copy_count_log_likelihoods: &[f64],
    mut current: &'a mut [f64],
    mut next: &'a mut [f64],
    log_allele_count_distribution: &mut [f64],
) {
    // **Both buffers, not just the one the first sample reads.** Each pass zeroes only the
    // window it is about to write — `[0, padding + live)` — and the *next* pass's `copies = 0`
    // tap reads `ploidy` entries beyond that, into the counts this sample has just made
    // reachable. Those entries are logically zero and must literally be zero, and after this
    // fill they stay so: no pass ever writes above its own window.
    //
    // **The alternative was measured and it is a silent wrong answer, not a crash.** Clearing
    // only `current` leaves `next` holding
    // [`UNWRITTEN_SCRATCH_VALUE`](crate::ng::calling::UNWRITTEN_SCRATCH_VALUE), which is
    // `NaN`; the `NaN` multiplies into the high counts from the second sample on, survives the
    // rescaling because `f64::max` returns the other operand, and is written out as `−∞`. The
    // count axis is then permanently truncated to `0..=ploidy` at every cohort size — at 63
    // diploid samples the site quality came out 46.3 Phred against production's 733.7 on the
    // same table.
    current.fill(0.0);
    next.fill(0.0);
    current[ploidy] = 1.0;
    let mut log_scale = 0.0;
    let mut log_at_zero = 0.0;

    for (sample, copy_counts) in copy_count_log_likelihoods
        .chunks_exact(ploidy + 1)
        .enumerate()
    {
        // The exact recurrence, independent of the linear bulk below.
        log_at_zero += copy_counts[0];

        // The row over its own maximum, so every tap lands in `(0, 1]`; the maximum folds
        // into the running scale.
        let largest = copy_counts
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let live = (sample + 1) * ploidy + 1;
        next[..ploidy + live].fill(0.0);
        for (copies, &log_weight) in copy_counts.iter().enumerate() {
            let weight = if log_weight > f64::NEG_INFINITY {
                (log_weight - largest).exp()
            } else {
                0.0
            };
            if weight == 0.0 {
                continue;
            }
            // Tap-major, so each pass is a contiguous multiply-add the compiler can
            // vectorise: `next[k] += current[k − copies] × weight`.
            let source = &current[(ploidy - copies)..(ploidy - copies) + live];
            let destination = &mut next[ploidy..ploidy + live];
            for (into, &from) in destination.iter_mut().zip(source) {
                *into += from * weight;
            }
        }

        log_scale += largest;
        let peak = next[ploidy..ploidy + live]
            .iter()
            .copied()
            .fold(0.0_f64, f64::max);
        if peak > 0.0 {
            let inverse = 1.0 / peak;
            for value in next[ploidy..ploidy + live].iter_mut() {
                *value *= inverse;
            }
            log_scale += peak.ln();
        }

        std::mem::swap(&mut current, &mut next);
    }

    for (count, slot) in log_allele_count_distribution.iter_mut().enumerate() {
        let value = current[ploidy + count];
        *slot = if value > 0.0 {
            value.ln() + log_scale
        } else {
            f64::NEG_INFINITY
        };
    }
    log_allele_count_distribution[0] = log_at_zero;
}

/// `ln Σ exp(v)` over a slice, in two passes: the largest entry, then the sum of the rest
/// relative to it. `−∞` where every entry is impossible.
fn log_sum_exp_over(values: &[f64]) -> f64 {
    let largest = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if largest == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let total: f64 = values.iter().map(|value| (value - largest).exp()).sum();
    largest + total.ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::genotype_prior::SeedRegime;
    use crate::ng::calling::{CallingScratch, CandidateAlleles, GenotypeTable};
    use crate::ng::locus_generation::LocusKind;
    use crate::ng::types::Ploidy;
    use std::sync::Arc;

    fn diploid() -> Ploidy {
        Ploidy::try_new(2).expect("a diploid")
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

    /// A neutral panel's fitted spectrum at one variant per kilobase — human diversity, and
    /// the middle column of `spec/calling_quality.md` §5.4's table.
    pub fn human_like_seed() -> SpectrumSeed {
        SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape)
    }

    /// Run the site quality over a cohort whose per-genotype log-likelihoods are `rows`.
    pub fn site_quality_of(rows: &[Vec<f64>], seed: SpectrumSeed) -> Phred {
        let (alleles, table) = generic_locus(1);
        let view = table.view();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(rows.len(), &alleles, &view);
        for (sample, row) in rows.iter().enumerate() {
            for (slot, &value) in scratch
                .sample_genotype_likelihoods_mut(sample)
                .iter_mut()
                .zip(row)
            {
                *slot = LogProb(value);
            }
        }
        score_uncorrected_site_quality(scratch.site_quality_buffers_mut(), &view, seed)
    }

    // -----------------------------------------------------------------------------------
    // The genotype quality
    // -----------------------------------------------------------------------------------

    /// **The whole formula on numbers a reader can follow.** A posterior of
    /// `[0.7, 0.2, 0.1]` puts 0.3 of the probability somewhere other than the winner, and
    /// `−10·log₁₀(0.3)` is 5.2288 Phred.
    #[test]
    fn one_samples_genotype_quality_matches_the_arithmetic_done_by_hand() {
        let (winner, quality) = score_best_genotype(&[0.7, 0.2, 0.1]);
        assert_eq!(winner, GenotypeIdx(0));
        assert!(
            (f64::from(quality.get()) - 5.228_787_452_803_376).abs() < 1e-4,
            "quality {} against -10 log10(0.3)",
            quality.get()
        );
    }

    /// **A sample whose reads leave one genotype standing does not produce an infinity.**
    /// `1 − 1.0` is zero and `log₁₀(0)` is `−∞`; the clamp to one unit in the last place
    /// below one is what makes the answer a number, and the cap is what makes it 99.
    ///
    /// **What the clamp is not doing, measured rather than assumed:** removing it
    /// (`let best_below_one = best;`) leaves this test and the whole module green. `-10 ·
    /// log₁₀(0)` is `+∞`, and `f32::clamp` maps `+∞` to the ceiling, so `Phred::try_new`
    /// never sees an infinity on *this* input. The clamp earns its place at
    /// `p_best` a hair below one, where the unclamped logarithm is finite but the clamped one
    /// differs — and against a future in which the cap is removed or raised. It is kept
    /// because it is production's and because `f32::clamp`'s treatment of an infinity is not
    /// a property this module should be resting on; **and it is currently untested**, which
    /// is worth knowing before anyone deletes it.
    #[test]
    fn a_certain_sample_is_capped_rather_than_infinite() {
        let (winner, quality) = score_best_genotype(&[0.0, 1.0, 0.0]);
        assert_eq!(winner, GenotypeIdx(1));
        assert_eq!(quality.get(), MAX_GENOTYPE_QUALITY);
    }

    /// **Two equally probable genotypes give the lower index**, because the fold keeps the
    /// first strict maximum. The genotype table's order is fixed, so a tie has one answer at
    /// any thread count — and this is the case where an implementation keeping the *last*
    /// maximum would differ.
    #[test]
    fn a_tie_goes_to_the_lower_genotype_index() {
        let (winner, quality) = score_best_genotype(&[0.5, 0.5]);
        assert_eq!(winner, GenotypeIdx(0));
        assert!(
            (f64::from(quality.get()) - 3.010_299_956_639_812).abs() < 1e-4,
            "half the probability elsewhere is -10 log10(0.5) = 3.01, got {}",
            quality.get()
        );
    }

    /// A row of no genotypes would return genotype 0 at a quality computed from `−∞` — a
    /// call at a locus with nothing to call.
    #[test]
    #[should_panic(expected = "a posterior row of none")]
    fn a_genotype_quality_over_no_genotypes_is_refused() {
        let _ = score_best_genotype(&[]);
    }

    /// A row that is entirely `NaN` never reaches a maximum at all — the fold ends holding
    /// the `−∞` it started from, which without a check clamps to a quality of zero on
    /// genotype 0 and reads as a confident reference call.
    #[test]
    #[should_panic(expected = "must total one")]
    fn a_posterior_row_that_is_entirely_nan_is_refused() {
        let _ = score_best_genotype(&[f64::NAN, f64::NAN, f64::NAN]);
    }

    /// **The shape the winner alone cannot see, and the reason the check is on the total.**
    /// One `NaN` beside two real probabilities loses every comparison, so the fold picks 0.7
    /// and returns an entirely ordinary-looking 5.229 Phred — a call made against a row one
    /// of whose genotypes has no probability at all. Summing catches it because `NaN`
    /// survives addition.
    #[test]
    #[should_panic(expected = "must total one")]
    fn a_single_nan_beside_real_probabilities_is_refused() {
        let _ = score_best_genotype(&[0.7, f64::NAN, 0.1]);
    }

    /// A row that does not sum to one is not a posterior. The E-step normalises it, so this
    /// is a statement about that normalisation rather than about a caller.
    #[test]
    #[should_panic(expected = "must total one")]
    fn a_posterior_row_that_does_not_sum_to_one_is_refused() {
        let _ = score_best_genotype(&[0.4, 1.5]);
    }

    // -----------------------------------------------------------------------------------
    // The site quality
    // -----------------------------------------------------------------------------------

    /// **At one sample the whole calculation has a closed form, and this is it, written
    /// independently.**
    ///
    /// One diploid sample gives a cohort allele count of 0, 1 or 2, so the fold is three
    /// terms and the posterior is `likelihood(k) × prior(k)` normalised over three entries.
    /// The oracle below builds the Beta-Binomial prior from the same seed but by its own
    /// arithmetic, so a defect in the collapse, the fold, the rescaling or the normalisation
    /// shows up as a disagreement.
    #[test]
    fn at_one_sample_the_site_quality_matches_the_closed_form() {
        let seed = human_like_seed();
        // Reads that favour the heterozygote by 3 nats over hom-ref and 5 over hom-alt.
        let row = vec![-3.0, 0.0, -5.0];
        let got = site_quality_of(std::slice::from_ref(&row), seed);

        // The oracle: collapse (at two alleles each genotype is already its own count),
        // multiply by the Beta-Binomial(2, α_alt, α_ref) prior, normalise, read entry 0.
        let alpha_alternative = seed.alpha_alt_total();
        let alpha_reference = seed.alpha_ref();
        let chromosomes = 2.0_f64;
        let log_beta_normaliser = lgamma(alpha_alternative) + lgamma(alpha_reference)
            - lgamma(alpha_alternative + alpha_reference);
        let unnormalised: Vec<f64> = (0..=2)
            .map(|count| {
                let count = f64::from(count);
                let log_ways = lgamma(chromosomes + 1.0)
                    - lgamma(count + 1.0)
                    - lgamma(chromosomes - count + 1.0);
                let log_split = lgamma(alpha_alternative + count)
                    + lgamma(alpha_reference + chromosomes - count);
                row[count as usize] + log_ways + log_split
                    - lgamma(alpha_alternative + alpha_reference + chromosomes)
                    - log_beta_normaliser
            })
            .collect();
        let largest = unnormalised
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let total: f64 = unnormalised.iter().map(|v| (v - largest).exp()).sum();
        let log_total = largest + total.ln();
        let want = -10.0 * (unnormalised[0] - log_total) / std::f64::consts::LN_10;

        assert!(
            (f64::from(got.get()) - want).abs() < 1e-3,
            "site quality {} against the closed form {want}",
            got.get()
        );
    }

    /// **The property the whole formula was chosen for: adding reference-looking samples
    /// does not pile up evidence for a variant.**
    ///
    /// The rejected alternative — multiply every sample's probability of being homozygous
    /// reference — is not a normalised posterior, so each additional reference-looking
    /// sample multiplies in another factor below one and its quality *climbs at a locus
    /// nobody carries* (`spec/calling_quality.md` §5.1). This runs both on the same fixture:
    /// one sample weakly favouring the heterozygote, and 0, 5 or 20 samples whose reads
    /// firmly say reference.
    ///
    /// **A locus nobody carries stays at nothing however large the cohort gets** — five
    /// hundred samples whose reads firmly say reference leave the quality at 1.3 in 100,000
    /// of a Phred.
    ///
    /// This is the property that matters operationally, and it is the one this module can
    /// demonstrate. **It is not the contrast `spec/calling_quality.md` §5.1 argues**, and the
    /// difference is recorded rather than papered over: §5.1 says the rejected
    /// `Π_s P(hom-ref)` formula *"grows with cohort size at a site nobody carries"* where the
    /// marginal stays bounded. Measured here, on the fixture below and on a thin one, **both
    /// formulas grow, and in the same proportion**: at 12 nats a cohort of 1 → 500 takes the
    /// rejected value from 0.0000267 to 0.0133 and the shipped one from 0.0000000267 to
    /// 0.0000134 — each 500 times its own base. At one nat the shipped quality grows
    /// *faster* than the rejected one (0.0019 → 831.88 against 1.77 → 885.11).
    ///
    /// **The growth has a legitimate cause the section does not mention**: in a cohort of 501
    /// thin samples, *nobody carries this* is a far stronger claim than in a cohort of one,
    /// and thin reads cannot support it. Whether §5.1's argument needs restating is the
    /// owner's, and it is raised in this step's report — the arithmetic is production's,
    /// unchanged, and it agrees with production to `1.2e-5` Phred.
    #[test]
    fn a_locus_nobody_carries_stays_at_nothing_however_many_samples_look_at_it() {
        let seed = human_like_seed();
        // Firmly reference: 12 nats against either genotype carrying an alternative copy.
        let reference_looking = vec![0.0, -12.0, -24.0];

        let alone =
            f64::from(site_quality_of(std::slice::from_ref(&reference_looking), seed).get());
        let crowd = f64::from(site_quality_of(&vec![reference_looking; 500], seed).get());
        assert!(
            crowd < 0.001,
            "five hundred samples that all say reference must not accumulate confidence in a \
             variant: the quality came out {crowd} Phred against {alone} at one sample"
        );
    }

    /// **The exact log-domain tracking of the entry at zero earns its place between 20 and
    /// 50 samples**, and this is the fixture that shows it.
    ///
    /// Fifty samples each favouring the heterozygote by 20 nats put `P(count = 0 | reads)` at
    /// about `e^−1000`, which no rescaled linear buffer holds. Measured, with the override
    /// (`log_allele_count_distribution[0] = log_at_zero`) deleted: 20 samples are unaffected
    /// at 1694.18, **50 samples go from 4295.97 to the 9999 ceiling**, and 80 from 6899.70 to
    /// 9999. So a cohort the size of the tomato panel is already past the point where
    /// removing it would report every strongly supported locus as maximally confident.
    ///
    /// **This test could not have been written before the fold was fixed.** With the fold
    /// reading uninitialised scratch, the count axis was truncated to `0..=ploidy` at every
    /// cohort size, the override never mattered, and the first version of this test recorded
    /// that as a finding about the device rather than about the defect.
    #[test]
    fn the_exact_zero_term_is_what_keeps_a_fifty_sample_cohort_off_the_ceiling() {
        let rows = vec![vec![-20.0, 0.0, -20.0]; 50];
        let quality = site_quality_of(&rows, human_like_seed());
        assert!(
            quality.get() < MAX_SITE_QUALITY,
            "hitting the ceiling means the entry at zero was read back off the underflowed \
             linear buffer rather than accumulated in logs: got {}",
            quality.get()
        );
        assert!(
            (4000.0..4600.0).contains(&quality.get()),
            "fifty samples at 20 nats each against a 1-in-1,000 spectrum came out at 4295.97 \
             Phred when this was written; {} is far enough away that the arithmetic changed \
             rather than drifted",
            quality.get()
        );
    }

    /// **A sample with no possible genotype is refused where its cause has a name.** Every
    /// entry of its collapsed row is `−∞`, the fold turns that into a `NaN`, and without this
    /// check the failure surfaces two steps later as a complaint about the normalisation.
    #[test]
    #[should_panic(expected = "candidate genotypes impossible")]
    fn a_sample_with_no_possible_genotype_is_refused() {
        let rows = vec![
            vec![-2.0, 0.0, -6.0],
            vec![f64::NEG_INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY],
        ];
        let _ = site_quality_of(&rows, human_like_seed());
    }

    /// A locus whose reads make the reference impossible gives a probability of zero, which
    /// the Phred scale has no number for. **Capping is this module's answer to
    /// [`Phred`]'s refusal**, and the ceiling is named once.
    #[test]
    fn a_reference_free_locus_is_capped_rather_than_refused() {
        let rows = vec![vec![f64::NEG_INFINITY, 0.0, -30.0]; 3];
        assert_eq!(
            site_quality_of(&rows, human_like_seed()).get(),
            MAX_SITE_QUALITY
        );
    }

    /// **A seed of zero alternative concentration does not make the prior a `NaN`.**
    /// [`SpectrumSeed`] admits exactly zero — a fully invariant cohort is a real answer —
    /// and `ln Γ(0)` is `+∞`, which would poison every term. The floor is what stops it.
    #[test]
    fn a_fully_invariant_seed_still_gives_a_quality() {
        let seed = SpectrumSeed::new(1.0, 0.0, SeedRegime::NeutralShape);
        let quality = site_quality_of(&[vec![-2.0, 0.0, -6.0]], seed);
        assert!(
            quality.get().is_finite(),
            "a seed with no alternative mass must still give a number, got {}",
            quality.get()
        );
    }

    /// **The collapse indexes the copy table by the allele count, not by the ploidy** — a
    /// distinction every other fixture in this file hides, because they are all diploid and
    /// biallelic, where the two numbers are both 2 and the two expressions coincide.
    ///
    /// Three alleles at ploidy 2 separates them: measured, reading the row at
    /// `genotype × ploidy` instead of `genotype × allele_count` moves this locus's quality
    /// from 4.843 to 0.274 Phred, and at ploidy 8 it indexes out of range.
    #[test]
    fn the_collapse_strides_by_the_allele_count_and_not_by_the_ploidy() {
        let (alleles, table) = generic_locus(2);
        let view = table.view();
        assert_eq!(view.allele_count(), 3);
        // The table's order at three alleles: 0/0, 0/1, 1/1, 0/2, 1/2, 2/2.
        let rows = [
            vec![-6.0, -2.0, 0.0, -6.0, -2.0, -6.0],
            vec![-6.0, -6.0, -6.0, -2.0, -2.0, 0.0],
            vec![0.0, -2.0, -6.0, -2.0, -6.0, -6.0],
        ];
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(rows.len(), &alleles, &view);
        for (sample, row) in rows.iter().enumerate() {
            for (slot, &value) in scratch
                .sample_genotype_likelihoods_mut(sample)
                .iter_mut()
                .zip(row)
            {
                *slot = LogProb(value);
            }
        }
        let quality = score_uncorrected_site_quality(
            scratch.site_quality_buffers_mut(),
            &view,
            human_like_seed(),
        );
        assert!(
            (f64::from(quality.get()) - 19.471_59).abs() < 1e-3,
            "a triallelic locus's quality came out at 19.472 Phred when this was written; \
             a collapse that strides by the ploidy gives 19.991. Got {}",
            quality.get()
        );
    }

    /// **Adding the same constant to every genotype's likelihood of every sample changes
    /// nothing**, because it multiplies every entry of the count distribution by the same
    /// factor and cancels in the normalisation.
    ///
    /// This is what exercises the running log scale: with every row peaking at exactly zero —
    /// which is how every other fixture here is written — the per-sample `log_scale +=
    /// largest` adds nothing on every pass, so deleting it changes no answer. Shifted rows
    /// make it carry a real number, and the invariance is the strongest statement available
    /// about the rescaling: the arithmetic must be blind to where the likelihoods sit.
    #[test]
    fn shifting_every_likelihood_by_a_constant_leaves_the_quality_alone() {
        let base = vec![
            vec![-2.0, 0.0, -6.0],
            vec![0.0, -1.0, -5.0],
            vec![-4.0, 0.0, -1.0],
        ];
        let unshifted = site_quality_of(&base, human_like_seed());
        for shift in [-37.0, 12.5] {
            let shifted: Vec<Vec<f64>> = base
                .iter()
                .map(|row| row.iter().map(|value| value + shift).collect())
                .collect();
            let got = site_quality_of(&shifted, human_like_seed());
            assert!(
                (f64::from(got.get()) - f64::from(unshifted.get())).abs() < 1e-4,
                "a shift of {shift} nats on every likelihood moved the quality from {} to \
                 {}, and it cancels in the normalisation so it must not",
                unshifted.get(),
                got.get()
            );
        }
    }

    /// **A quality above the ceiling but finite is capped too**, not only an infinite one.
    /// Four hundred samples each favouring the heterozygote by 20 nats reach about 34,690
    /// Phred, which is a real number the `QUAL` column has no use for.
    #[test]
    fn a_finite_quality_above_the_ceiling_is_capped() {
        let rows = vec![vec![-20.0, 0.0, -20.0]; 400];
        assert_eq!(
            site_quality_of(&rows, human_like_seed()).get(),
            MAX_SITE_QUALITY
        );
    }

    /// **The prior is the run's own spectrum, and swapping it moves the answer in the
    /// direction §5.4 predicts.** On a panel ten times as polymorphic, a locus being variant
    /// is more likely before any read is looked at, so the same reads earn a higher quality.
    #[test]
    fn a_more_polymorphic_seed_earns_the_same_reads_a_higher_quality() {
        let rows = vec![vec![-2.0, 0.0, -6.0], vec![0.0, -6.0, -12.0]];
        let human = site_quality_of(
            &rows,
            SpectrumSeed::new(1.0, 1e-3, SeedRegime::NeutralShape),
        );
        let diverse = site_quality_of(
            &rows,
            SpectrumSeed::new(1.0, 1e-2, SeedRegime::NeutralShape),
        );
        assert!(
            diverse.get() > human.get(),
            "a tenfold more polymorphic seed must not lower the quality of the same reads: \
             {} against {}",
            diverse.get(),
            human.get()
        );
    }

    /// A cohort of no samples would fold nothing and read a quality off an empty axis.
    #[test]
    #[should_panic(expected = "a site quality over none")]
    fn a_site_quality_over_no_samples_is_refused() {
        let (alleles, table) = generic_locus(1);
        let view = table.view();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(1, &alleles, &view);
        let mut buffers = scratch.site_quality_buffers_mut();
        buffers.sample_count = 0;
        let _ = score_uncorrected_site_quality(buffers, &view, human_like_seed());
    }

    /// The likelihood table must be samples × genotypes. `chunks_exact` and the row slicing
    /// below it would otherwise read another sample's row, or panic on an index that names
    /// neither the cohort nor the table.
    #[test]
    #[should_panic(expected = "the genotype likelihoods are samples × genotypes")]
    fn a_likelihood_table_of_the_wrong_shape_is_refused() {
        let (alleles, table) = generic_locus(1);
        let view = table.view();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &view);
        let mut buffers = scratch.site_quality_buffers_mut();
        buffers.genotype_likelihoods = &buffers.genotype_likelihoods[..3];
        let _ = score_uncorrected_site_quality(buffers, &view, human_like_seed());
    }

    /// The count axis is padded by ploidy on the left so the fold's inner loop needs no
    /// bounds test. An unpadded one would slice out of range mid-fold, naming an offset
    /// rather than the buffer that is the wrong size.
    #[test]
    #[should_panic(expected = "padded by ploidy on the left")]
    fn an_unpadded_count_axis_is_refused() {
        let (alleles, table) = generic_locus(1);
        let view = table.view();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(2, &alleles, &view);
        let mut spare = vec![0.0; 5];
        let mut buffers = scratch.site_quality_buffers_mut();
        buffers.allele_count_distribution = &mut spare;
        let _ = score_uncorrected_site_quality(buffers, &view, human_like_seed());
    }

    /// The collapsed rows are one per sample, `ploidy + 1` wide. A short table would leave
    /// the last samples' rows unwritten — and `chunks_exact` walks only whole rows, so the
    /// fold would sum a cohort smaller than the one that was scored, silently.
    #[test]
    #[should_panic(expected = "the copy-count log-likelihoods are samples")]
    fn a_copy_count_table_of_the_wrong_shape_is_refused() {
        let (alleles, table) = generic_locus(1);
        let view = table.view();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(3, &alleles, &view);
        let mut spare = vec![0.0; 6];
        let mut buffers = scratch.site_quality_buffers_mut();
        buffers.copy_count_log_likelihoods = &mut spare;
        let _ = score_uncorrected_site_quality(buffers, &view, human_like_seed());
    }

    /// The log-domain result runs from a count of zero to `ploidy × samples` inclusive. One
    /// entry short and the prior is applied to a shorter axis than the fold filled, so the
    /// normalisation divides by a total that is missing its own tail.
    #[test]
    #[should_panic(expected = "the count axis runs from 0 to ploidy")]
    fn a_log_axis_of_the_wrong_length_is_refused() {
        let (alleles, table) = generic_locus(1);
        let view = table.view();
        let mut scratch = CallingScratch::<()>::default();
        scratch.prepare_for_locus(3, &alleles, &view);
        let mut spare = vec![0.0; 4];
        let mut buffers = scratch.site_quality_buffers_mut();
        buffers.log_allele_count_distribution = &mut spare;
        let _ = score_uncorrected_site_quality(buffers, &view, human_like_seed());
    }

    /// **A row that sums to one and still is not a distribution.** `[1.5, -0.5]` totals
    /// exactly one, so the total check passes it; the winning probability is 1.5, which is
    /// not a probability, and this is the only check that sees it.
    #[test]
    #[should_panic(expected = "must be finite and within [0, 1]")]
    fn a_winning_probability_above_one_is_refused_even_where_the_row_totals_one() {
        let _ = score_best_genotype(&[1.5, -0.5]);
    }
}
