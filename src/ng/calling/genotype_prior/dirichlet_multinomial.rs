//! The log-prior row: a handful of allele copies drawn from a population whose composition
//! is not known. That is the Dirichlet-multinomial — the locus's allele frequencies drawn
//! from a Dirichlet and averaged out, with this sample's copies drawn from what is left.
//!
//! **This file holds the random-mating half.** Every genotype gets it. The other half, the
//! branch for two copies that are the same one counted twice, is the inbreeding mixture that
//! plan step B2 wraps around it (`doc/devel/ng/spec/calling_priors.md` §3.2).

use crate::genetics::{PROBABILITY_FLOOR, lgamma};
use crate::ng::calling::genotype_prior::{GenotypePriorModel, PriorRow};
use crate::ng::types::{InbreedingF, LogProb};

/// Fill the row with each genotype's **random-mating** log-prior — what the genotype would be
/// worth if the sample's copies were independent draws from the population.
///
/// For a genotype carrying `k_a` copies of allele `a`, summing to the ploidy `m`:
///
/// ```text
/// log P_random(g) = log C(m; k) + Σ_a [ lgamma(α_a + k_a) − lgamma(α_a) ]
/// ```
///
/// `log C(m; k)` is how many orderings of the genome's copies spell the genotype — `ln 2` for
/// a diploid heterozygote, `ln 1 = 0` for a homozygote.
///
/// ## Up to a shared additive constant, and why that is not a shortcut
///
/// The exact Dirichlet-multinomial also carries `lgamma(Σα + m) − lgamma(Σα)`, which is the
/// same for every genotype. It is **omitted**: it cancels the moment the loop rescales the
/// row, so carrying it would cost two `lgamma` calls per sample per pass to move every entry
/// by one number. What comes out is therefore not a normalised log-probability — an entry may
/// be positive, and the row sums to nothing in particular (spec §3.1).
///
/// ## Where it comes from, and which of production's two spellings this is
///
/// A port of `dirichlet_multinomial_log_priors` (`src/genetics.rs`), the shared primitive the
/// plan and the architecture name as this step's source. **One thing changed: it fills the
/// caller's row, and the per-allele `lgamma(α_a)` baseline goes in the caller's scratch, where
/// production allocated a `Vec` for each.** Nothing may allocate inside the per-sample loop
/// (spec §8).
///
/// **Production has that shared primitive and a second, private copy of the same mathematics,
/// and they do not agree to the last bit.** The shared one has exactly one shipping caller,
/// the STR cohort's EM (`src/ssr/cohort/em.rs`); the SNP/indel engine runs
/// `fill_log_indep_per_g_from` (`src/var_calling/posterior_engine.rs`) instead, which already
/// takes a caller's `out` and a caller's `lgamma_alpha` — the same shape this port needed,
/// arrived at independently and for the same stated reason. It **associates differently**,
/// summing the per-allele terms and adding the multinomial coefficient last where this folds
/// from the coefficient. Measured over a grid of 492 genotype values, **112 of them differ, by
/// at most one unit in the last place**.
///
/// So the parity test below pins this port to the *shared* primitive, and anyone building a
/// differential against production's **SNP/indel** path should expect an ulp of disagreement
/// from a function this file does not copy. It cannot move a genotype — a log-prior ulp against
/// read likelihoods that differ by whole nats — but the GIAB 83.6% → 94.6% measurement of spec
/// §2.2 was taken on that other path, so the two are worth telling apart before anyone
/// reconciles them.
///
/// ## Preconditions
///
/// All structural: held by [`PriorRow::new`], in release. The value precondition — every
/// concentration entry finite and at least the alternative floor — is [`super::Concentration`]'s,
/// and it is **debug-only**, which is where production holds the same check. In release a
/// concentration entry of zero gives `lgamma(0) = +∞` and a row entry of `−∞`, and a `NaN`
/// entry gives a `NaN` row; neither is detected here, and both mean a seed builder handed over
/// something it should have floored.
pub fn fill_random_mating_log_priors(row: &mut PriorRow<'_>) {
    let concentration = row.concentration().get();
    let allele_count = concentration.len();
    let genotype_allele_counts = row.genotype_allele_counts();
    let log_multinomial_coeffs = row.log_multinomial_coeffs();
    let (lgamma_concentration, out) = row.scratch_and_out();

    // lgamma(α_a) once per allele — the baseline every genotype's term subtracts. This is the
    // whole reason the seam carries a per-allele scratch: recomputing it at each (genotype,
    // allele) pair would nearly double the `lgamma` count.
    for (slot, &alpha) in lgamma_concentration.iter_mut().zip(concentration) {
        *slot = lgamma(alpha);
    }

    for ((slot, copies_of), &log_coeff) in out
        .iter_mut()
        .zip(genotype_allele_counts.chunks_exact(allele_count))
        .zip(log_multinomial_coeffs)
    {
        *slot = LogProb(one_genotypes_log_prior(
            copies_of,
            log_coeff,
            concentration,
            lgamma_concentration,
        ));
    }
}

/// `log C(m; k) + Σ_a [ lgamma(α_a + k_a) − lgamma(α_a) ]` for one genotype, where
/// `copies_of[a]` is `k_a` and `lgamma_concentration[a]` is the `lgamma(α_a)` the caller
/// computed once.
///
/// **The order of the additions is part of the contract, not an implementation detail.** The
/// fold starts at `log_coeff` and takes each carried allele's two terms in allele order, which
/// is what the primitive this ports does. Summing the alleles first and adding the coefficient
/// last — the tidier spelling, and the one production's *other* copy uses — moves a row by one
/// unit in the last place and fails the bit-parity test.
///
/// **An allele the genotype carries no copy of is skipped, and that is not only a saving.** Its
/// term `lgamma(α_a + 0) − lgamma(α_a)` is zero in exact arithmetic, so the loop does one
/// `lgamma` per allele the genotype actually carries rather than one per allele in the table.
/// But adding a zero term is not free in floating point: the fold associates as
/// `(acc + lgamma(α_a + k_a)) − lgamma(α_a)`, so a skipped allele's two large, nearly equal
/// logarithms would enter and leave the accumulator with a rounding in between. Measured —
/// removing the branch moves a diploid biallelic hom-ref row from `0.6931471805599453` to
/// `0.6931471805599454`.
#[inline]
fn one_genotypes_log_prior(
    copies_of: &[u32],
    log_coeff: f64,
    concentration: &[f64],
    lgamma_concentration: &[f64],
) -> f64 {
    copies_of
        .iter()
        .zip(concentration)
        .zip(lgamma_concentration)
        .fold(log_coeff, |acc, ((&copies, &alpha), &lgamma_alpha)| {
            if copies == 0 {
                acc
            } else {
                acc + lgamma(alpha + f64::from(copies)) - lgamma_alpha
            }
        })
}

/// The default answer at step 8's seam: §3.1's Dirichlet-multinomial under §3.2's two-branch
/// inbreeding mixture.
///
/// **Two branches, not a correction term.** The inbreeding coefficient `F` is the chance that a
/// sample's two copies of a locus are *identical by descent* — one ancestral copy counted twice —
/// rather than two independent draws from the population. So the prior is a mixture over which of
/// those happened:
///
/// ```text
/// with probability F        the copies are one copy counted twice, so the genotype is
///                           homozygous for allele a with probability α_a / Σα
/// with probability 1 − F    the copies are independent draws, and the genotype is worth
///                           what fill_random_mating_log_priors gave it
/// ```
///
/// Only the second branch exists for a genotype that is not homozygous, so its row entry is just
/// `log(1 − F)` plus the random-mating term; a homozygous one takes the `logsumexp` of both.
///
/// **Why the mixture rather than the textbook Wright formulas.** `P(het) = 2pq(1 − F)` and
/// `p² + Fpq` are biallelic and diploid. This says the same thing at two alleles and keeps saying
/// something at three, at four, and at a tetraploid — which is why the Wright formulas are a test
/// oracle here and not a code path (spec §3.2).
///
/// ## The two branches have to be on the same scale, and one correction puts them there
///
/// [`fill_random_mating_log_priors`] returns log-priors **up to a shared additive constant** —
/// it drops `lgamma(Σα + m) − lgamma(Σα)` because it is the same for every genotype and cancels
/// when the loop rescales the row. That is true of a row on its own and **false the moment a
/// second branch is mixed into it**: the identical-by-descent term `α_a / Σα` is a true
/// probability, so mixing the two directly inflates the random-mating branch by `Σα(Σα + 1)` at
/// diploidy and the inbreeding coefficient does a fraction of the work it should.
///
/// **Measured, biallelic diploid, at tomato1's fitted diversity of 6 in 10,000.** Read the
/// heterozygote-to-homozygous-alternative prior ratio, which is what the coefficient is for. An
/// outbred sample sits at 2:1; the more inbred the sample, the further below that it should fall.
///
/// ```text
///                     outbred    what the model says    what the uncorrected
///                     (F = 0)    at F = 0.8             mixture gives at F = 0.8
///   1 sample             2.00           0.222                    0.400
///   50 samples         188.7            0.493                  181.8
///   1,000 samples     1818              0.499                 1816.5
/// ```
///
/// At **one sample** the uncorrected mixture makes a heterozygote 1.80 times as likely as it
/// should be (1.90 at `F = 0.9`) — the ratio still travels 90% of the way from the outbred value
/// to the right one. **At cohort scale it barely moves at all**: 50 samples leave it 3.6% of the
/// way and 1,000 samples 0.09%, so an inbreeding coefficient of 0.8 buys almost nothing. The
/// reason is that the concentration this function is handed is the leave-one-out one (spec §6), so
/// `Σα` grows with the cohort and the inflation factor `Σα(Σα + 1)` grows with its square.
///
/// With the correction the row matches Wright's formulas to **1.9e-5** in the concentrated limit,
/// inside the `1e-4` the test below allows.
///
/// **Production has the defect and it is live rather than latent.** Its engine mixes the same
/// two scales (`posterior_engine.rs`), and its own *default* inbreeding coefficient is `0`, where
/// the branch short-circuits away and nothing shows — but the pipeline also hands the engine the
/// per-sample coefficients the diversity estimator **fitted**, as overrides, and those are not
/// zero on an inbred panel (`var_calling/pipeline.rs`, `with_fixation_index_overrides`). ng
/// corrects it deliberately (owner, 2026-08-22); this is the one place the port departs from what
/// it was ported from, and spec §3.1's "the constant cancels" needs the qualification that it
/// cancels in a row and not in a mixture.
///
/// ## Where this term stops being computable
///
/// `lgamma(Σα + m) − lgamma(Σα)` subtracts two nearly equal numbers of order `Σα·ln Σα`, so it
/// loses precision as `Σα` grows. Measured as the row's departure from one unit of probability,
/// biallelic, ploidy 2 and 8, `F` up to 0.95: **1.5e-11 at `Σα = 7.2e3`** — about 3,600 diploid
/// samples, past the top of the committed cohort range — rising to 9.1e-11 at 1.2e5, **1.1e-9 at
/// 1.2e6** and 2.8e-7 at 1.2e8. Nothing in the caller's range is affected; the figure at 1.2e6 is
/// worth knowing because [`the_concentrated_limit_matches_the_wright_formulas`] drives `Σα` there
/// deliberately, which is why the normalisation identity is not also run at that total.
///
/// **The correction is added to the identical-by-descent branch rather than subtracted from the
/// other**, which is the same mixture up to the shared constant a row is allowed to carry, and
/// leaves an outbred sample's row bit-identical to the primitive's.
///
/// A port of the mixture in `src/var_calling/posterior_engine.rs`, in its two-branch form.
///
/// **Derives `Debug` like every other public type in this module**, so a `&dyn GenotypePriorModel`
/// can be printed: the seam exists to compare two priors, and a result that cannot name which one
/// produced it is not auditable after the fact (spec §2.2).
#[derive(Copy, Clone, PartialEq, Debug)]
pub struct MarginalizedDirichletPrior;

impl GenotypePriorModel for MarginalizedDirichletPrior {
    fn name(&self) -> &'static str {
        "marginalized-dirichlet"
    }

    fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, inbreeding: InbreedingF) {
        fill_inbreeding_mixture_log_priors(row, inbreeding.get());
    }
}

/// The mixture on a bare coefficient rather than on the newtype.
///
/// **It exists so `F = 1` can be tested**, which is the mathematical edge of the model and not a
/// case the caller is meant to meet: production's estimator clamps at 0.99, and ng's newtype is to
/// be tightened to `[0, 1)` by the prerequisites plan. **That tightening has not happened** —
/// [`InbreedingF::try_new(1.0)`](crate::ng::types::InbreedingF::try_new) returns `Ok` today — so
/// `F = 1` is reachable through the trait as well, and both spellings are pinned by tests. The
/// limit is worth pinning either way: at `F = 1` every heterozygote becomes impossible, and the two
/// homozygotes must stand in the ratio `α_ref : α_alt` (spec §7, §12 test 3).
///
/// **This is not a test-only path**, whatever its reason for existing: the trait implementation
/// above routes every caller through it, which is why the coefficient is checked here rather than
/// left to the newtype. The check is `debug_assert!`, which is where this module puts every *value*
/// check (`Concentration::new`); the structural checks that guard silent truncation are the ones
/// held in release.
fn fill_inbreeding_mixture_log_priors(row: &mut PriorRow<'_>, inbreeding: f64) {
    debug_assert!(
        (0.0..=1.0).contains(&inbreeding),
        "the inbreeding coefficient must be a fraction in [0, 1]; got {inbreeding}. A value \
         outside it — or a NaN — makes the whole row NaN rather than failing, and a negative one \
         poisons only the homozygotes, which normalises to a plausible wrong answer"
    );

    fill_random_mating_log_priors(row);

    let concentration = row.concentration().get();
    let homozygous_allele_for = row.homozygous_allele_for();
    let concentration_total = concentration.iter().sum::<f64>();
    // ln(Σα), so the identical-by-descent branch's α_a / Σα costs one `ln` per homozygous
    // genotype rather than a division and a logarithm.
    let log_concentration_total = concentration_total.ln();
    // The genotype-independent term the primitive drops, `lgamma(Σα + m) − lgamma(Σα)`, put
    // back — on this branch rather than taken off the other, so that at `F = 0` the row is the
    // primitive's bit for bit and nothing downstream of an outbred sample shifts. See the
    // function's doc for why dropping it is safe in a row and not in a mixture.
    let shared_normalising_term =
        lgamma(concentration_total + f64::from(row.ploidy())) - lgamma(concentration_total);
    // The two mixture weights, and they are floored differently on purpose.
    //
    // `1 − F` is floored because it is the one that can reach a *row entry*: at `F = 1` every
    // heterozygote's only branch is this one, so an unfloored `ln(0)` would write `−∞` into the
    // output. Spec §8 and arch §1.1 both say a genotype the prior rules out carries a finite, very
    // negative log-prior rather than `−∞`, and the comparator arriving at plan step F1 is ported
    // from `wright_genotype_log_priors`, which floors internally — so this is also what keeps the
    // two implementations behind the seam on one convention.
    //
    // `F` is *not* floored, because its `ln(0)` at `F = 0` never reaches a row entry: it makes the
    // identical-by-descent branch impossible, and `log_sum_exp_2` returns the other branch exactly.
    // Flooring it would replace that short-circuit with two `exp` and a `ln` per homozygous
    // genotype on every outbred sample — the ordinary case — to move nothing.
    let log_weight_identical_by_descent = inbreeding.ln();
    let log_weight_independent_draws = (1.0 - inbreeding).max(PROBABILITY_FLOOR).ln();

    let (_, out) = row.scratch_and_out();
    // `zip` would truncate to the shorter of the two, which is the silent failure this module's
    // own doc names as the worst available here — `PriorRow::new` holds their equality in release
    // so it cannot happen.
    for (slot, homozygous) in out.iter_mut().zip(homozygous_allele_for) {
        let independent_draws = log_weight_independent_draws + slot.get();
        *slot = LogProb(match homozygous {
            // The homozygous test is this lookup and nothing else — no comparison of the copy
            // counts anywhere in this function. That is what gives the above-diploidy question one
            // place to change (spec §3.3).
            Some(allele) => {
                let identical_by_descent = log_weight_identical_by_descent
                    + concentration[usize::from(allele.0)].ln()
                    - log_concentration_total
                    + shared_normalising_term;
                log_sum_exp_2(independent_draws, identical_by_descent)
            }
            None => independent_draws,
        });
    }
}

/// `ln(e^a + e^b)`, shifted by the larger so neither exponential overflows.
///
/// **The two `−∞` short-circuits are a saving and a floor — not a correctness guard on any row
/// this caller can produce.** An earlier version of this comment claimed otherwise and was wrong;
/// it is recorded here because the claim reads as load-bearing and is cheap to believe. Measured:
/// delete either short-circuit, or both, and every entry of every fixture in this module is
/// bit-identical. The general path already returns the finite argument exactly when the other is
/// `−∞`, because `a.max(b)` picks the finite one and `ln(1 + 0)` is `0`.
///
/// What they do buy is two things:
///
/// - **A saving on the ordinary case, and it is the second guard that provides it.** At `F = 0` the
///   identical-by-descent branch is `−∞` on every homozygous genotype, so the guard on `b` fires
///   and skips two `exp` and a `ln` per homozygote. **This is why the mixture floors `1 − F` and
///   not `F`**: flooring both would put every outbred sample down the general path to move nothing.
/// - **A pair of `−∞` arguments returns `−∞` rather than `NaN`.** `−∞ − −∞` is the only route to a
///   `NaN` here. Since the mixture floors `1 − F`, argument `a` is now never `−∞`, so the first
///   guard cannot fire from any coefficient and the pair needs a concentration entry of zero, which
///   [`Concentration::new`](super::Concentration) refuses. The guard stays because this is a shared
///   helper and the next caller may not floor.
///
/// Ported from the same engine as the mixture, whose own comment gives the saving as the reason.
///
/// **Visible to the folder because both implementations behind the seam use it**: the comparator
/// in [`hardy_weinberg`](super::hardy_weinberg) mixes the same two branches over a different
/// random-mating term, and a second spelling of this would let the two disagree about `−∞` without
/// anything saying so.
#[inline]
pub(super) fn log_sum_exp_2(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let larger = a.max(b);
    larger + ((a - larger).exp() + (b - larger).exp()).ln()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ng::calling::GenotypeTable;
    use crate::ng::calling::genotype_prior::Concentration;
    use crate::ng::types::Ploidy;

    /// `ln` of the rising factorial `Π_{j=0}^{k−1} (a + j)` — the closed form of
    /// `lgamma(a + k) − lgamma(a)` for an integer `k`, **computed with no `lgamma` at all**.
    ///
    /// Carried across from production's own test module (`src/genetics.rs`). It is an
    /// independent implementation rather than a table of golden values, which is what lets it
    /// keep checking after a constant moves or the concentration is fitted differently.
    fn pochhammer_ln(alpha: f64, copies: u32) -> f64 {
        (0..copies).map(|j| (alpha + f64::from(j)).ln()).sum()
    }

    /// One genotype's log-prior the second way: `log_coeff + Σ_a pochhammer_ln(α_a, k_a)`.
    fn rising_factorial_log_prior(copies_of: &[u32], log_coeff: f64, concentration: &[f64]) -> f64 {
        log_coeff
            + copies_of
                .iter()
                .zip(concentration)
                .map(|(&copies, &alpha)| pochhammer_ln(alpha, copies))
                .sum::<f64>()
    }

    /// A concentration of the shape the loop hands over: the reference's entry, and
    /// `alternative_total` **shared unevenly across the alternatives so that they sum to it** —
    /// uneven because equal entries would let a fold that paired alleles with the wrong
    /// concentration still come out right.
    ///
    /// `reference` is a parameter because the array this primitive is handed is not the seed but
    /// the leave-one-out concentration (spec §6): the seed plus the cohort's expected allele
    /// copies. At one sample that leaves the reference near 1; at a thousand diploid samples it
    /// is near 2,000, and the parity test below runs both.
    fn concentration_of(allele_count: usize, reference: f64, alternative_total: f64) -> Vec<f64> {
        let mut values = vec![reference; allele_count];
        let weight_total: f64 = (1..allele_count).map(|index| index as f64).sum();
        for (index, slot) in values.iter_mut().enumerate().skip(1) {
            *slot = alternative_total * index as f64 / weight_total;
        }
        values
    }

    /// Run the primitive over one shape and hand back the row and the concentration it used.
    fn row_for(
        copies: u8,
        allele_count: usize,
        reference: f64,
        alternative_total: f64,
    ) -> (Vec<LogProb>, Vec<f64>) {
        let table = GenotypeTable::build(Ploidy::try_new(copies).unwrap(), allele_count);
        let view = table.view();
        let concentration = concentration_of(allele_count, reference, alternative_total);
        let mut scratch = vec![0.0; allele_count];
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut row = PriorRow::new(
            Concentration::new(&concentration),
            view.genotype_allele_counts(),
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &mut scratch,
            &mut out,
        );
        fill_random_mating_log_priors(&mut row);
        (out, concentration)
    }

    /// **The independent oracle** (spec §12 test 4): every log-prior matches a rising-factorial
    /// computation using no `lgamma` at all, over shapes from a haploid monomorphic locus to an
    /// octoploid one, at four diversities apiece.
    ///
    /// **The `1e-12` tolerance is a fact about this grid rather than about the primitive**, and
    /// saying so is the point. It is production's own figure for the same comparison, and it
    /// holds here because every concentration on this grid has its reference entry at 1. The
    /// true disagreement grows with that entry — measured at 7.2e-14 when it is 201, and
    /// **2.1e-12 when it is 2,001**, which is a thousand diploid samples' worth of leave-one-out
    /// counts. So this test stays at the small end deliberately, and the **bit-parity** test
    /// below covers the large one, where there is no tolerance to argue about.
    ///
    /// Nor is this the test that catches a re-associated fold: that moves a row by one unit in
    /// the last place and sails through here.
    #[test]
    fn every_log_prior_matches_the_rising_factorial_oracle() {
        let mut rows_checked = 0_usize;
        for (copies, allele_count) in [
            (1_u8, 1_usize),
            (1, 4),
            (2, 1),
            (2, 2),
            (2, 3),
            (2, 6),
            (3, 4),
            (4, 3),
            (6, 2),
            (8, 2),
        ] {
            for alternative_total in [1e-4, 1e-3, 1e-2, 0.5] {
                let (row, concentration) = row_for(copies, allele_count, 1.0, alternative_total);
                let table = GenotypeTable::build(Ploidy::try_new(copies).unwrap(), allele_count);
                let view = table.view();
                for (genotype, (copies_of, &log_coeff)) in view
                    .genotype_allele_counts()
                    .chunks_exact(allele_count)
                    .zip(view.log_multinomial_coeffs())
                    .enumerate()
                {
                    let want = rising_factorial_log_prior(copies_of, log_coeff, &concentration);
                    assert!(
                        (row[genotype].get() - want).abs() < 1e-12,
                        "ploidy {copies}, {allele_count} alleles, alternative total \
                         {alternative_total}, genotype {genotype} {copies_of:?}: got {} want \
                         {want}",
                        row[genotype].get()
                    );
                    rows_checked += 1;
                }
            }
        }
        // 87 genotypes over the ten shapes — 1, 4, 1, 3, 6, 21, 20, 15, 7 and 9, each
        // `C(alleles + copies − 1, copies)` — at four diversities apiece. Asserted so that a
        // shape silently dropped from the grid fails here rather than quietly narrowing what
        // the oracle covers.
        assert_eq!(
            rows_checked, 348,
            "the grid should cover 348 genotype rows; it covered {rows_checked}"
        );
    }

    /// **The port agrees with what it was ported from, bit for bit.**
    ///
    /// The oracle above proves the mathematics; this proves the arithmetic is *performed the
    /// same way*, which is a different claim and the one a port owes. Bit equality rather than a
    /// tolerance, because the only route to it is the same operations in the same order — which
    /// a re-associated fold quietly breaks by one unit in the last place.
    ///
    /// **The grid runs the reference entry from 1 to 6,001, and that is the load-bearing part.**
    /// What this primitive is handed is the leave-one-out concentration, the run's seed plus the
    /// cohort's expected allele copies (spec §6), so the reference entry is near 1 only at one
    /// sample and near 2,000 at a thousand diploid ones. With the grid pinned at 1 — as it was
    /// when this test was first written — three separate wrong implementations passed the whole
    /// module: one clamping every entry at 1, one hard-coding the reference's entry, and one
    /// truncating the fold at six alleles. They move a row by 9.9 nats at a hundred samples and
    /// get 57 of 78 genotypes wrong on a twelve-allele locus.
    ///
    /// Reading production here is an oracle and never a dependency: nothing ng ships imports
    /// from `src/genetics.rs` beyond `lgamma` and the alternative-concentration floor.
    #[test]
    fn the_port_matches_production_bit_for_bit() {
        for (copies, allele_count) in [
            (1_u8, 3_usize),
            (2, 1),
            (2, 2),
            (2, 4),
            (2, 12),
            (3, 3),
            (4, 2),
            (8, 2),
        ] {
            for reference in [1.0, 201.0, 2001.0, 6001.0] {
                for alternative_total in [1e-3, 1e-2, 0.25, 40.0] {
                    let (row, concentration) =
                        row_for(copies, allele_count, reference, alternative_total);
                    let table =
                        GenotypeTable::build(Ploidy::try_new(copies).unwrap(), allele_count);
                    let view = table.view();
                    let production = crate::genetics::dirichlet_multinomial_log_priors(
                        view.genotype_allele_counts(),
                        view.log_multinomial_coeffs(),
                        allele_count,
                        &concentration,
                    );
                    assert_eq!(production.len(), row.len());
                    for (genotype, (&ours, &theirs)) in row.iter().zip(&production).enumerate() {
                        assert_eq!(
                            ours.get().to_bits(),
                            theirs.to_bits(),
                            "ploidy {copies}, {allele_count} alleles, reference {reference}, \
                             alternative total {alternative_total}, genotype {genotype}: ours {} \
                             theirs {theirs}",
                            ours.get()
                        );
                    }
                }
            }
        }
    }

    /// **An allele the genotype carries no copy of is skipped, and the skip is bit-exact.**
    ///
    /// Proved here without reading production, which is what makes it worth keeping: the parity
    /// test above is otherwise the only thing that catches the skip going away, and it stops
    /// being able to the day `src/genetics.rs` moves.
    ///
    /// The hom-reference genotype at a locus with one allele and the hom-reference genotype at a
    /// locus with four both carry two copies of allele 0 and no copy of anything else, so the
    /// two must agree to the last bit. Add the skipped alleles back as zero terms and they do
    /// not: each pair of large, nearly equal logarithms enters and leaves the accumulator with a
    /// rounding in between — the one unit in the last place this file's doc measures.
    #[test]
    fn an_allele_a_genotype_does_not_carry_cannot_move_its_prior() {
        let (one_allele, _) = row_for(2, 1, 1.0, 0.0);
        for alternative_total in [1e-4, 1e-3, 1e-2, 0.5] {
            let (four_alleles, concentration) = row_for(2, 4, 1.0, alternative_total);
            assert_eq!(
                four_alleles[0].get().to_bits(),
                one_allele[0].get().to_bits(),
                "hom-reference moved when three alleles it carries no copy of were added \
                 (concentration {concentration:?}): {} against {}",
                four_alleles[0].get(),
                one_allele[0].get()
            );
            // Genotype 1/1 does carry allele 1, so it must differ from hom-reference —
            // otherwise the assertion above would pass on a function that ignored the
            // concentration entirely.
            assert_ne!(
                four_alleles[2].get(),
                four_alleles[0].get(),
                "1/1 carries allele 1 and must not land on hom-reference's value"
            );
        }
    }

    /// **The 2:1 tripwire** (spec §12 test 1). At a biallelic diploid locus with no inbreeding
    /// the heterozygote is twice as likely as the homozygous-alternative genotype — exactly
    /// `2·α_ref : (1 + α_alt)`, which is what the neutral spectrum written as a Dirichlet gives,
    /// and which reaches 2:1 as `α_alt → 0`.
    ///
    /// **This is spec §2.3's trap seen in a test.** Production's plug-in path regularised its
    /// frequency estimate with `α_ref = 10`, and the *marginalized* prior over that same
    /// Dirichlet gives 20:1 — the same wrong answer computed more expensively. (The 22:1 the
    /// spec also records is a different quantity: the plug-in path's own Hardy–Weinberg ratio at
    /// the frequency that regularisation implies. The two are easy to conflate and are not the
    /// same number.)
    ///
    /// **What this guards, and what it does not.** The reference entry it uses is written in
    /// this test module, so nothing here stops a *seed builder* from choosing 10 — that is step
    /// D3's to guard. What it pins is that this function turns the pair it is given into the
    /// ratio the mathematics requires, and the last assertion shows it moving when the pair
    /// does.
    #[test]
    fn the_heterozygote_is_twice_the_homozygous_alternative_at_every_realistic_diversity() {
        for alternative_total in [1e-4, 6e-4, 1e-3, 1e-2] {
            let (row, _) = row_for(2, 2, 1.0, alternative_total);
            let ratio = (row[1].get() - row[2].get()).exp();
            let want = 2.0 / (1.0 + alternative_total);
            assert!(
                (ratio / want - 1.0).abs() < 1e-12,
                "at α_alt {alternative_total} the het:hom-alt ratio was {ratio}, not {want}"
            );
        }

        // The widest diversity above is 1 in 100, where the exact ratio is 2/1.01 = 1.9802, so
        // 2% is the band the whole loop fits inside and 1% is not. Asserted at that one point
        // rather than inside the loop, so adding a wilder diversity cannot silently widen it.
        let (widest, _) = row_for(2, 2, 1.0, 1e-2);
        let widest_ratio = (widest[1].get() - widest[2].get()).exp();
        assert!(
            (widest_ratio - 2.0).abs() < 0.02,
            "at 1 in 100 the ratio should still be within 2% of 2:1, got {widest_ratio}"
        );

        // And the same function at a reference of 1.5 is not near 2:1 at all, which is what
        // makes the assertions above a check on the input rather than an identity.
        let (raised, _) = row_for(2, 2, 1.5, 1e-3);
        let raised_ratio = (raised[1].get() - raised[2].get()).exp();
        assert!(
            (raised_ratio - 2.997).abs() < 1e-3,
            "at a reference of 1.5 the ratio should be about 2.997, got {raised_ratio}"
        );
    }

    /// **The invariant mass follows the diversity** (spec §12 test 2): the homozygous-reference
    /// genotype takes about `1 − 3θ/2` of the row, so raising the diversity moves mass onto the
    /// genotypes that carry a variant.
    ///
    /// Checked against the closed form rather than a remembered number, at a tolerance that
    /// scales with `θ²` — the neglected term — so the assertion tightens as the approximation
    /// improves rather than passing on slack. It is genuinely tight: the true error is about
    /// `1.75θ²` against a `3θ²` budget, so each of the three uses a little under 60% of it.
    ///
    /// **There is no monotonicity assertion, because it could not fail.** The `1 − 1.5θ ± 3θ²`
    /// windows at these three θ do not overlap, so the tolerance check already forces the
    /// ordering; an assertion that the weights fall would have read like a second check and been
    /// none.
    #[test]
    fn the_homozygous_reference_weight_follows_the_diversity() {
        for theta in [1e-4, 1e-3, 1e-2] {
            let (row, _) = row_for(2, 2, 1.0, theta);
            let largest = row
                .iter()
                .map(|p| p.get())
                .fold(f64::NEG_INFINITY, f64::max);
            let total: f64 = row.iter().map(|p| (p.get() - largest).exp()).sum();
            let hom_reference_weight = (row[0].get() - largest).exp() / total;

            assert!(
                (hom_reference_weight - (1.0 - 1.5 * theta)).abs() < 3.0 * theta * theta,
                "at θ {theta} the hom-ref weight was {hom_reference_weight}, not about {}",
                1.0 - 1.5 * theta
            );
        }
    }

    /// Run the mixture over one shape and hand back the row it wrote.
    fn mixed_row_for(
        copies: u8,
        allele_count: usize,
        reference: f64,
        alternative_total: f64,
        inbreeding: f64,
    ) -> Vec<LogProb> {
        let table = GenotypeTable::build(Ploidy::try_new(copies).unwrap(), allele_count);
        let view = table.view();
        let concentration = concentration_of(allele_count, reference, alternative_total);
        let mut scratch = vec![0.0; allele_count];
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut row = PriorRow::new(
            Concentration::new(&concentration),
            view.genotype_allele_counts(),
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &mut scratch,
            &mut out,
        );
        fill_inbreeding_mixture_log_priors(&mut row, inbreeding);
        out
    }

    /// **At no inbreeding the mixture is the identity on the random-mating row**, bit for bit.
    ///
    /// It has to be: `log(1 − 0)` is exactly zero, and the identical-by-descent branch enters at
    /// `log 0 = −∞`, which the `logsumexp` short-circuit returns the other argument for. So an
    /// outbred sample pays nothing for the branch existing, and — the reason this is asserted
    /// rather than argued — everything the previous tests pin about the primitive continues to
    /// hold through the seam, including the 2:1 ratio and the invariant mass.
    ///
    /// Bit equality rather than a tolerance, because anything else would mean the mixture is
    /// computing rather than short-circuiting, and the `+ 0.0` alone would flip a `−0.0`.
    #[test]
    fn at_no_inbreeding_the_mixture_leaves_the_random_mating_row_untouched() {
        for (copies, allele_count) in [(1_u8, 2_usize), (2, 2), (2, 3), (2, 6), (4, 3), (8, 2)] {
            for alternative_total in [1e-4, 1e-3, 1e-2, 0.5] {
                let (random_mating, _) = row_for(copies, allele_count, 1.0, alternative_total);
                let mixed = mixed_row_for(copies, allele_count, 1.0, alternative_total, 0.0);
                for (genotype, (mixed, random)) in mixed.iter().zip(&random_mating).enumerate() {
                    assert_eq!(
                        mixed.get().to_bits(),
                        random.get().to_bits(),
                        "ploidy {copies}, {allele_count} alleles, alternative total \
                         {alternative_total}, genotype {genotype}: mixed {} against random-mating \
                         {}",
                        mixed.get(),
                        random.get()
                    );
                }
            }
        }
    }

    /// **The mixed row is a probability distribution once its shared constant is removed** — at
    /// every ploidy, every allele count and every inbreeding coefficient.
    ///
    /// This is the identity the two branches being on one scale *is*. The mixture says the
    /// genotype came either from two independent draws or from one copy counted twice, so the
    /// weights have to add to one:
    ///
    /// ```text
    /// Σ_g P(g) = (1 − F)·Σ_g P_random(g) + F·Σ_a (α_a / Σα) = (1 − F) + F = 1
    /// ```
    ///
    /// The row carries a shared `lgamma(Σα + m) − lgamma(Σα)` on top, because the correction is
    /// added to the identical-by-descent branch rather than taken off the other; subtract it and
    /// the row must sum to exactly one.
    ///
    /// **It is the only test here that is not diploid**, and that is deliberate: the Wright oracle
    /// is biallelic-diploid by construction, so without this one a correction that hard-coded the
    /// ploidy at 2 passes the whole module — measured, it did. It is also what would catch the
    /// correction being dropped again.
    ///
    /// **The reference entry runs to 6,001, and that is the other load-bearing part.** What the
    /// mixture is handed is the leave-one-out concentration (spec §6) — near 1 at one sample and
    /// near 2,000 at a thousand diploid ones — and at 1.0 alone a `Σα` that hard-codes the
    /// reference entry is indistinguishable from the true one. Measured: such an implementation
    /// leaves this identity at 4.0e-15 with the reference pinned at 1.0, and misses by 0.95 once
    /// it reaches 6,001. This is the same narrowing that B1's parity grid was widened to fix.
    ///
    /// Worst error over the whole grid on the correct code is **1.5e-11**, against the `1e-9`
    /// budget. The budget is not widened for the large totals because the term's own precision is
    /// what sets it: past `Σα` of about a million the identity would need a looser one (see
    /// [`fill_inbreeding_mixture_log_priors`]), and no cohort in range reaches that.
    #[test]
    fn the_mixed_row_is_a_probability_distribution_once_its_shared_constant_is_removed() {
        for (copies, allele_count) in [
            (1_u8, 1_usize),
            (1, 3),
            (2, 1),
            (2, 2),
            (2, 4),
            (3, 3),
            (4, 2),
            (4, 4),
            (8, 2),
        ] {
            for inbreeding in [0.0, 0.25, 0.8, 0.95] {
                for alternative_total in [1e-3, 1e-2, 0.5] {
                    for reference in [1.0, 201.0, 2001.0, 6001.0] {
                        let concentration =
                            concentration_of(allele_count, reference, alternative_total);
                        let total: f64 = concentration.iter().sum();
                        let shared_constant = crate::genetics::lgamma(total + f64::from(copies))
                            - crate::genetics::lgamma(total);
                        let row = mixed_row_for(
                            copies,
                            allele_count,
                            reference,
                            alternative_total,
                            inbreeding,
                        );
                        let mass: f64 = row.iter().map(|p| (p.get() - shared_constant).exp()).sum();
                        assert!(
                            (mass - 1.0).abs() < 1e-9,
                            "ploidy {copies}, {allele_count} alleles, F {inbreeding}, reference \
                             {reference}, alternative total {alternative_total}: the row carries \
                             {mass} rather than one unit of probability"
                        );
                    }
                }
            }
        }
    }

    /// **The Wright oracle** (spec §3.2, §12): a second independent check of the mixture, from
    /// formulas that share no code with it.
    ///
    /// The two meet in a limit. A Dirichlet with a large total is concentrated on one frequency
    /// `p = α_alt / Σα`, so the Dirichlet-multinomial collapses to a plain binomial draw at that
    /// `p` — and the two-branch mixture then says exactly what Wright's biallelic diploid formulas
    /// say: `2pq(1 − F)` for the heterozygote, `q² + Fpq` and `p² + Fpq` for the homozygotes. So
    /// the concentration is scaled by a million **at a fixed ratio**, which moves the total without
    /// moving `p`.
    ///
    /// Compared as *differences between genotypes*, because this module's rows are log-priors up to
    /// a shared additive constant while `wright_genotype_log_priors` returns normalised ones. Run
    /// at `F = 0` and `F = 0.5`, the two the spec names, and at three frequencies.
    #[test]
    fn the_concentrated_limit_matches_the_wright_formulas() {
        for frequency in [0.05, 0.2, 0.5] {
            for inbreeding in [0.0, 0.5] {
                // Σα = 1e6 at the wanted p, so the Dirichlet is effectively a point mass there.
                let total = 1e6;
                let alternative = total * frequency;
                let reference = total - alternative;
                let row = mixed_row_for(2, 2, reference, alternative, inbreeding);

                let (hom_reference, heterozygote, hom_alternative) =
                    crate::genetics::wright_genotype_log_priors(frequency, inbreeding);

                // Row order is the VCF one: 0/0, 0/1, 1/1.
                let ours_het_over_hom_ref = row[1].get() - row[0].get();
                let ours_hom_alt_over_hom_ref = row[2].get() - row[0].get();
                let wright_het_over_hom_ref = heterozygote - hom_reference;
                let wright_hom_alt_over_hom_ref = hom_alternative - hom_reference;

                assert!(
                    (ours_het_over_hom_ref - wright_het_over_hom_ref).abs() < 1e-4,
                    "p {frequency}, F {inbreeding}: het over hom-ref was \
                     {ours_het_over_hom_ref}, Wright says {wright_het_over_hom_ref}"
                );
                assert!(
                    (ours_hom_alt_over_hom_ref - wright_hom_alt_over_hom_ref).abs() < 1e-4,
                    "p {frequency}, F {inbreeding}: hom-alt over hom-ref was \
                     {ours_hom_alt_over_hom_ref}, Wright says {wright_hom_alt_over_hom_ref}"
                );
            }
        }
    }

    /// **The `1e-4` tolerance above is the limit's own error, and it closes as the limit is
    /// approached** — which is what makes the oracle a check rather than a coincidence.
    ///
    /// The Dirichlet-multinomial reaches Wright's formulas only as `Σα → ∞`; at a finite total the
    /// gap is of order `1/Σα`. Measured at four totals a decade apart: each tenfold rise in the
    /// concentration shrinks the disagreement by about tenfold, so a tolerance that passed at
    /// `Σα = 1e6` for the wrong reason would not narrow with it.
    ///
    /// **The first total carries an assertion too**, which it did not when the accumulator started
    /// at `f64::INFINITY` — three of the four comparisons were real and the first was
    /// `gap < INFINITY`. It is seeded instead with a measured bound: the four gaps are 1.10e-2,
    /// 1.11e-3, 1.11e-4 and 1.11e-5, and 1e-1 puts the first threshold at 2e-2, comfortably above
    /// the 1.10e-2 the correct code produces and far below the **7.98e-1** the uncorrected mixture
    /// produces there — so the seed is what makes the dropped correction fail on the first
    /// iteration rather than the second.
    #[test]
    fn the_wright_agreement_closes_as_the_concentration_grows() {
        let frequency = 0.2;
        let inbreeding = 0.5;
        let (hom_reference, heterozygote, _) =
            crate::genetics::wright_genotype_log_priors(frequency, inbreeding);
        let wright_het_over_hom_ref = heterozygote - hom_reference;

        // Seeded with a measured bound rather than `f64::INFINITY`, so the first total carries an
        // assertion too. See this test's doc for the four gaps and why 1e-1 is the seed.
        let mut previous_gap = 1e-1;
        let mut budget_is_a_seed = true;
        for total in [1e2, 1e3, 1e4, 1e5] {
            let alternative = total * frequency;
            let row = mixed_row_for(2, 2, total - alternative, alternative, inbreeding);
            let gap = ((row[1].get() - row[0].get()) - wright_het_over_hom_ref).abs();
            let comparand = if budget_is_a_seed {
                "the seeded bound"
            } else {
                "the gap at a tenth the concentration"
            };
            assert!(
                gap < previous_gap / 5.0,
                "at Σα {total} the gap to Wright was {gap}, not a fifth of {previous_gap} \
                 ({comparand}) — the limit is not being approached"
            );
            previous_gap = gap;
            budget_is_a_seed = false;
        }
    }

    /// **The full-inbreeding limit** (spec §12 test 3). At `F = 1` every heterozygote is impossible
    /// and the two homozygotes stand in the ratio `α_ref : α_alt` — the identical-by-descent branch
    /// alone, which draws one allele and counts it twice.
    ///
    /// **The estimator never delivers this**, so it tests the mathematics at its edge rather than a
    /// case the caller meets: production clamps its inbreeding estimate at 0.99, and ng's newtype is
    /// *to be* tightened to `[0, 1)`. **It has not been** — `InbreedingF::try_new(1.0)` returns `Ok`
    /// today, which is why the sibling test drives the same limit through the seam as well. Once the
    /// tightening lands, the largest coefficient the type can carry is `1 − 2⁻⁵³`, where `1 − F` is
    /// about `1.1e-16` and the floor below never bites; this path keeps the limit reachable, because
    /// the bare-coefficient function admits `1` by design.
    ///
    /// **The heterozygote lands at the probability floor, not at `−∞`** (spec §8, arch §1.1, owner
    /// 2026-08-22). `log(1 − 1)` is `ln 0`, and the mixture floors the weight at
    /// [`PROBABILITY_FLOOR`] first, so the entry is about `−691` plus the genotype's own
    /// random-mating value — finite, and 300 orders of magnitude below anything a read can move.
    /// Asserted as a bound rather than an exact number so the concentration may vary.
    #[test]
    fn at_full_inbreeding_only_the_homozygotes_survive_and_they_stand_at_the_concentration_ratio() {
        for alternative_total in [1e-3, 1e-2, 0.25] {
            let row = mixed_row_for(2, 2, 1.0, alternative_total, 1.0);
            assert!(
                row[1].get().is_finite() && row[1].get() < -600.0,
                "the heterozygote must be impossible at F = 1 but finite, got {}",
                row[1].get()
            );
            let ratio = (row[0].get() - row[2].get()).exp();
            assert!(
                (ratio - 1.0 / alternative_total).abs() / (1.0 / alternative_total) < 1e-12,
                "hom-ref : hom-alt was {ratio}, not the concentration ratio {}",
                1.0 / alternative_total
            );
        }
    }

    /// **The homozygous branch fires on the table's lookup and on nothing else** (spec §3.3).
    ///
    /// This is the property the whole above-diploidy question rests on: what "homozygous" should
    /// mean when four copies can be two identical-by-descent and two not is deferred to a spec of
    /// its own, and its cost of change is one function only while nothing recomputes the test from
    /// the copy counts.
    ///
    /// Checked by lying to it. The row is built over a real diploid biallelic table but handed a
    /// lookup that says no genotype is homozygous; both homozygotes must then come back with the
    /// independent-draws branch alone, exactly as the heterozygote does. An implementation that
    /// compared copy counts against the ploidy — the obvious inline spelling — would ignore the lie
    /// and fail here.
    #[test]
    fn the_homozygous_branch_reads_the_tables_lookup_rather_than_the_copy_counts() {
        let table = GenotypeTable::build(Ploidy::try_new(2).unwrap(), 2);
        let view = table.view();
        let concentration = [1.0, 1e-2];
        let inbreeding = 0.6;

        let truthful = mixed_row_for(2, 2, 1.0, 1e-2, inbreeding);

        let mut scratch = [0.0; 2];
        let mut lied_to = vec![LogProb(f64::NAN); view.genotype_count()];
        let no_genotype_is_homozygous = vec![None; view.genotype_count()];
        let mut row = PriorRow::new(
            Concentration::new(&concentration),
            view.genotype_allele_counts(),
            view.log_multinomial_coeffs(),
            &no_genotype_is_homozygous,
            &mut scratch,
            &mut lied_to,
        );
        fill_inbreeding_mixture_log_priors(&mut row, inbreeding);

        let (random_mating, _) = row_for(2, 2, 1.0, 1e-2);
        let log_outbreeding = (1.0 - inbreeding).ln();
        for (genotype, entry) in lied_to.iter().enumerate() {
            assert_eq!(
                entry.get().to_bits(),
                (log_outbreeding + random_mating[genotype].get()).to_bits(),
                "with the lookup saying nothing is homozygous, genotype {genotype} should carry \
                 the independent-draws branch alone"
            );
        }
        // And the truthful lookup does move the two homozygotes, so the assertion above is not
        // passing because the inbreeding branch does nothing at this concentration.
        assert_ne!(truthful[0].get(), lied_to[0].get());
        assert_ne!(truthful[2].get(), lied_to[2].get());
        assert_eq!(truthful[1].get().to_bits(), lied_to[1].get().to_bits());
    }

    /// **Raising the inbreeding coefficient moves mass from the heterozygote onto the
    /// homozygotes**, which is the whole of what the coefficient is for.
    ///
    /// Stated as a ratio rather than as two weights so it does not depend on the row's shared
    /// additive constant, and checked against the closed form the mixture implies at a
    /// concentration low enough that the independent-draws branch dominates: the heterozygote
    /// carries a factor `(1 − F)` that neither homozygote does, so the ratio falls like `(1 − F)`.
    #[test]
    fn raising_the_inbreeding_coefficient_moves_mass_onto_the_homozygotes() {
        let alternative_total = 1e-2;
        let mut previous = f64::INFINITY;
        for inbreeding in [0.0, 0.2, 0.5, 0.8, 0.95] {
            let row = mixed_row_for(2, 2, 1.0, alternative_total, inbreeding);
            let het_over_hom_alt = (row[1].get() - row[2].get()).exp();
            assert!(
                het_over_hom_alt < previous,
                "at F {inbreeding} the het:hom-alt ratio was {het_over_hom_alt}, not below the \
                 {previous} at the lower coefficient"
            );
            previous = het_over_hom_alt;
        }
        // At F = 0.95 the heterozygote has lost twenty-fold against random mating, and the
        // homozygous-alternative genotype has gained the identical-by-descent branch, so the ratio
        // is far under the 2:1 an outbred sample sees.
        let outbred = mixed_row_for(2, 2, 1.0, alternative_total, 0.0);
        assert!(
            (outbred[1].get() - outbred[2].get()).exp() > 1.9,
            "an outbred sample should still see about 2:1"
        );
        assert!(
            previous < 0.2,
            "at F = 0.95 the het:hom-alt ratio should be far under 2:1, got {previous}"
        );
    }

    /// Run the mixture **through the seam**, exactly as every caller outside this file will.
    fn seam_row_for(
        copies: u8,
        allele_count: usize,
        reference: f64,
        alternative_total: f64,
        inbreeding: InbreedingF,
    ) -> Vec<LogProb> {
        let table = GenotypeTable::build(Ploidy::try_new(copies).unwrap(), allele_count);
        let view = table.view();
        let concentration = concentration_of(allele_count, reference, alternative_total);
        let mut scratch = vec![0.0; allele_count];
        let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
        let mut row = PriorRow::new(
            Concentration::new(&concentration),
            view.genotype_allele_counts(),
            view.log_multinomial_coeffs(),
            view.homozygous_alleles(),
            &mut scratch,
            &mut out,
        );
        MarginalizedDirichletPrior.fill_genotype_log_priors(&mut row, inbreeding);
        out
    }

    /// **The trait implementation is the only thing outside this file can reach, and every other
    /// test here drives the private function instead** — so this is what ties the two together.
    ///
    /// The three lines of the implementation are where the coefficient could be dropped, clamped or
    /// inverted, and none of it is visible to a test that calls the mixture directly. Measured: an
    /// implementation that ignored the coefficient and called the random-mating primitive left the
    /// whole module green while moving the het:hom-alt ratio at `F = 0.95` from 0.051 to 1.980 —
    /// 39-fold, about 16 on the Phred scale.
    ///
    /// Bit equality against the bare-coefficient spelling, so the two cannot drift apart, plus one
    /// assertion that the coefficient reached the mixture at all.
    #[test]
    fn the_seam_and_the_bare_coefficient_agree_and_both_carry_the_inbreeding_coefficient() {
        for (copies, allele_count) in [(1_u8, 2_usize), (2, 2), (2, 4), (4, 3)] {
            for inbreeding in [0.0, 0.2, 0.8, 0.95] {
                let through_the_seam = seam_row_for(
                    copies,
                    allele_count,
                    1.0,
                    1e-2,
                    InbreedingF::try_new(inbreeding).unwrap(),
                );
                let directly = mixed_row_for(copies, allele_count, 1.0, 1e-2, inbreeding);
                for (genotype, (seam, direct)) in through_the_seam.iter().zip(&directly).enumerate()
                {
                    assert_eq!(
                        seam.get().to_bits(),
                        direct.get().to_bits(),
                        "ploidy {copies}, {allele_count} alleles, F {inbreeding}, genotype \
                         {genotype}: the seam wrote {} and the bare coefficient {}",
                        seam.get(),
                        direct.get()
                    );
                }
            }
        }

        // And the coefficient is not merely passed but used: at a biallelic diploid locus an
        // implementation that dropped it would leave the outbred 2:1 ratio in place.
        let inbred = seam_row_for(2, 2, 1.0, 1e-2, InbreedingF::try_new(0.95).unwrap());
        let ratio = (inbred[1].get() - inbred[2].get()).exp();
        assert!(
            ratio < 0.1,
            "at F = 0.95 the het:hom-alt ratio through the seam should be far under 2:1, got \
             {ratio}"
        );
    }

    /// **The full-inbreeding limit on the path a caller actually takes.**
    ///
    /// [`at_full_inbreeding_only_the_homozygotes_survive_and_they_stand_at_the_concentration_ratio`]
    /// drives the bare coefficient because [`InbreedingF`] is *to be* tightened to `[0, 1)`. It has
    /// not been — measured, `InbreedingF::try_new(1.0)` returns `Ok` on this commit — so the limit
    /// is reachable through the seam today and is pinned there too.
    ///
    /// **When the tightening lands this constructor stops returning `Ok`**, and that is the signal
    /// to move this test to the largest coefficient the newtype then accepts rather than to delete
    /// it.
    #[test]
    fn the_seam_rules_out_heterozygotes_at_the_greatest_coefficient_the_newtype_accepts() {
        let one = InbreedingF::try_new(1.0)
            .expect("InbreedingF still accepts 1.0; tighten this test when it stops");
        for alternative_total in [1e-3, 1e-2, 0.25] {
            let row = seam_row_for(2, 2, 1.0, alternative_total, one);
            assert!(
                row[1].get().is_finite() && row[1].get() < -600.0,
                "the heterozygote must be impossible at F = 1 but finite, got {}",
                row[1].get()
            );
            assert!(
                row[0].get().is_finite() && row[2].get().is_finite(),
                "both homozygotes must stay finite: {:?}",
                row.iter().map(|p| p.get()).collect::<Vec<_>>()
            );
            // Every entry finite is the property the floor exists for: a row the normaliser can
            // subtract a maximum from without producing a NaN, whatever the coefficient.
            assert!(
                row.iter().all(|p| p.get().is_finite()),
                "no entry may be −∞ once the weight is floored: {:?}",
                row.iter().map(|p| p.get()).collect::<Vec<_>>()
            );
            let ratio = (row[0].get() - row[2].get()).exp();
            assert!(
                (ratio - 1.0 / alternative_total).abs() / (1.0 / alternative_total) < 1e-12,
                "hom-ref : hom-alt was {ratio}, not the concentration ratio {}",
                1.0 / alternative_total
            );
        }
    }

    /// **The `−∞` guards, pinned directly**, because no row fixture can fail on them.
    ///
    /// Measured: delete either guard, or both, and every entry of every row in this module is
    /// bit-identical — the unguarded path computes `a + ln(1 + 0)`, which is `a` exactly. What the
    /// guards genuinely decide is the both-`−∞` pair, which without them is `−∞ − −∞` and so a
    /// `NaN`. That pair means both branches of the mixture are impossible, and the answer to it is
    /// an impossible genotype, not a missing number.
    #[test]
    fn log_sum_exp_2_returns_the_finite_argument_when_the_other_is_impossible() {
        assert_eq!(log_sum_exp_2(f64::NEG_INFINITY, -3.5), -3.5);
        assert_eq!(log_sum_exp_2(-3.5, f64::NEG_INFINITY), -3.5);
        assert_eq!(
            log_sum_exp_2(f64::NEG_INFINITY, f64::NEG_INFINITY),
            f64::NEG_INFINITY,
            "two impossible branches make an impossible genotype, not a NaN"
        );
        // The ordinary case is still the log-sum-exp: ln(e^0 + e^0) = ln 2.
        assert!((log_sum_exp_2(0.0, 0.0) - std::f64::consts::LN_2).abs() < 1e-15);
        // And the shift by the larger argument is what keeps it finite when one side is far away:
        // without it `exp(-800)` underflows to zero and `exp(800)` overflows to infinity.
        assert!((log_sum_exp_2(-800.0, -800.0) - (-800.0 + std::f64::consts::LN_2)).abs() < 1e-12);
    }

    /// **At haploidy the coefficient can do nothing, and every entry goes through the mixture.**
    ///
    /// One copy is one draw, so "the copies are one ancestral copy counted twice" and "the copies
    /// are independent draws" are the same statement, and the mixture must return the random-mating
    /// row whatever `F` is. Every haploid genotype is homozygous, which makes this the only test
    /// where the identical-by-descent branch fires on *every* genotype — a stronger check on the
    /// two branches being on one scale than the row adding to one, because they have to coincide
    /// entry by entry rather than in total.
    ///
    /// A tolerance rather than bit equality: measured, `F = 0.95` moves an entry by about 1.4e-17,
    /// the rounding of adding two logarithms that cancel.
    #[test]
    fn a_haploid_row_is_unmoved_by_the_inbreeding_coefficient() {
        for allele_count in [1_usize, 2, 3, 6] {
            for alternative_total in [1e-4, 1e-2, 0.5] {
                let (random_mating, _) = row_for(1, allele_count, 1.0, alternative_total);
                for inbreeding in [0.0, 0.25, 0.8, 0.95] {
                    let mixed = mixed_row_for(1, allele_count, 1.0, alternative_total, inbreeding);
                    for (genotype, (mixed, random)) in mixed.iter().zip(&random_mating).enumerate()
                    {
                        assert!(
                            (mixed.get() - random.get()).abs() < 1e-12,
                            "{allele_count} alleles, alternative total {alternative_total}, F \
                             {inbreeding}, genotype {genotype}: haploid mixture {} against \
                             random-mating {}",
                            mixed.get(),
                            random.get()
                        );
                    }
                }
            }
        }
    }

    /// **A locus with one allele has one genotype and `α_0 / Σα = 1`, so both branches say the same
    /// thing** — the row is the primitive's whatever `F` is.
    ///
    /// The degenerate case of the correction being exactly right, at every ploidy, and a shape the
    /// other mixture tests reach only at `F = 0`. A tolerance rather than bit equality: measured,
    /// `F = 0.95` moves it by one unit in the last place.
    #[test]
    fn a_monomorphic_locus_is_unmoved_by_the_inbreeding_coefficient() {
        for copies in [1_u8, 2, 4, 8] {
            let (random_mating, _) = row_for(copies, 1, 1.0, 0.0);
            assert_eq!(random_mating.len(), 1);
            for inbreeding in [0.0, 0.5, 0.95, 1.0] {
                let mixed = mixed_row_for(copies, 1, 1.0, 0.0, inbreeding);
                assert!(
                    (mixed[0].get() - random_mating[0].get()).abs() < 1e-12,
                    "ploidy {copies}, F {inbreeding}: a monomorphic locus moved from {} to {}",
                    random_mating[0].get(),
                    mixed[0].get()
                );
            }
        }
    }

    /// **`ploidy()` reads one genotype's copy counts rather than a stored number**, so this pins
    /// both halves of that: that the slice it reads is exactly one genotype wide at every shape —
    /// including the one-allele table, where over-reaching would run off the end of a one-row
    /// table — and that the premise it rests on holds, every genotype summing to the same count.
    ///
    /// The mixture adds `lgamma(Σα + m) − lgamma(Σα)` to one branch only, so a wrong `m` is not a
    /// shared constant: it re-weights the mixture.
    #[test]
    fn ploidy_returns_the_copy_count_every_genotype_sums_to() {
        for copies in [1_u8, 2, 3, 4, 6, 8] {
            for allele_count in [1_usize, 2, 3, 6] {
                let table = GenotypeTable::build(Ploidy::try_new(copies).unwrap(), allele_count);
                let view = table.view();
                let concentration = concentration_of(allele_count, 1.0, 1e-2);
                let mut scratch = vec![0.0; allele_count];
                let mut out = vec![LogProb(f64::NAN); view.genotype_count()];
                let row = PriorRow::new(
                    Concentration::new(&concentration),
                    view.genotype_allele_counts(),
                    view.log_multinomial_coeffs(),
                    view.homozygous_alleles(),
                    &mut scratch,
                    &mut out,
                );
                assert_eq!(
                    row.ploidy(),
                    u32::from(copies),
                    "ploidy {copies}, {allele_count} alleles"
                );
                for (genotype, copies_of) in view
                    .genotype_allele_counts()
                    .chunks_exact(allele_count)
                    .enumerate()
                {
                    assert_eq!(
                        copies_of.iter().sum::<u32>(),
                        u32::from(copies),
                        "ploidy {copies}, {allele_count} alleles, genotype {genotype}"
                    );
                }
            }
        }
    }
}
