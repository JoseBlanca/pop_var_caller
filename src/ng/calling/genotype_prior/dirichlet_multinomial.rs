//! The log-prior row: a handful of allele copies drawn from a population whose composition
//! is not known. That is the Dirichlet-multinomial — the locus's allele frequencies drawn
//! from a Dirichlet and averaged out, with this sample's copies drawn from what is left.
//!
//! **This file holds the random-mating half.** Every genotype gets it. The other half, the
//! branch for two copies that are the same one counted twice, is the inbreeding mixture that
//! plan step B2 wraps around it (`doc/devel/ng/spec/calling_priors.md` §3.2).

use crate::genetics::lgamma;
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
/// **Measured at the concentration this caller ships.** One sample at tomato1's fitted diversity
/// of 6 in 10,000, biallelic diploid: without the correction the heterozygote-to-hom-alt prior
/// ratio is 0.400 at `F = 0.8` where the model says 0.222, and 0.200 at `F = 0.9` against 0.105
/// — a heterozygote made about 1.8 times as likely as it should be, in the direction this caller
/// is already weakest. With it, the row matches the Wright formulas to one part in a million in
/// the concentrated limit, which is the test below.
///
/// **Production has the defect and it is live rather than latent.** Its engine mixes the same
/// two scales (`posterior_engine.rs`), and its own default inbreeding coefficient is `0`, where
/// the branch short-circuits away and nothing shows — but its pipeline passes the cohort's
/// *fitted* coefficient, so any run on an inbred panel meets it. ng corrects it deliberately
/// (owner, 2026-08-22); this is the one place the port departs from what it was ported from, and
/// spec §3.1's "the constant cancels" needs the qualification that it cancels in a row and not
/// in a mixture.
///
/// **The correction is added to the identical-by-descent branch rather than subtracted from the
/// other**, which is the same mixture up to the shared constant a row is allowed to carry, and
/// leaves an outbred sample's row bit-identical to the primitive's.
///
/// A port of the mixture in `src/var_calling/posterior_engine.rs`, in its two-branch form.
pub struct MarginalizedDirichletPrior;

impl GenotypePriorModel for MarginalizedDirichletPrior {
    fn fill_genotype_log_priors(&self, row: &mut PriorRow<'_>, inbreeding: InbreedingF) {
        fill_marginalized_log_priors(row, inbreeding.get());
    }
}

/// The mixture on a bare coefficient rather than on the newtype.
///
/// **It exists so `F = 1` can be tested**, which is the mathematical edge of the model and not a
/// case the caller meets: production's estimator clamps at 0.99, and ng's newtype is to be tightened
/// to `[0, 1)` by the prerequisites plan, so nothing will be able to hand the trait a 1. The limit
/// is still worth pinning — at `F = 1` every heterozygote becomes impossible, and the two
/// homozygotes must stand in the ratio `α_ref : α_alt` (spec §7, §12 test 3).
fn fill_marginalized_log_priors(row: &mut PriorRow<'_>, inbreeding: f64) {
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
    let scale_of_the_random_branch =
        lgamma(concentration_total + f64::from(row.ploidy())) - lgamma(concentration_total);
    // `ln(0)` is `−∞` and that is the wanted value in both directions: at `F = 0` the
    // identical-by-descent branch is impossible, at `F = 1` the independent-draws branch is.
    let log_inbreeding = inbreeding.ln();
    let log_outbreeding = (1.0 - inbreeding).ln();

    let (_, out) = row.scratch_and_out();
    for (slot, homozygous) in out.iter_mut().zip(homozygous_allele_for) {
        let independent_draws = log_outbreeding + slot.get();
        *slot = LogProb(match homozygous {
            // The homozygous test is this lookup and nothing else — no comparison of the copy
            // counts anywhere in this function. That is what gives the above-diploidy question one
            // place to change (spec §3.3).
            Some(allele) => {
                let identical_by_descent = log_inbreeding
                    + concentration[usize::from(allele.0)].ln()
                    - log_concentration_total
                    + scale_of_the_random_branch;
                log_sum_exp_2(independent_draws, identical_by_descent)
            }
            None => independent_draws,
        });
    }
}

/// `ln(e^a + e^b)`, shifted by the larger so neither exponential overflows.
///
/// **Both short-circuits are load-bearing rather than defensive.** At `F = 0` the
/// identical-by-descent branch is `−∞` on every homozygous genotype and at `F = 1` the
/// independent-draws branch is `−∞` on every genotype, so these are the ordinary cases and not
/// edge ones — and without the first, `−∞ − −∞` would make every row `NaN`. Ported from the same
/// engine as the mixture.
#[inline]
fn log_sum_exp_2(a: f64, b: f64) -> f64 {
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
        fill_marginalized_log_priors(&mut row, inbreeding);
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
                    let concentration = concentration_of(allele_count, 1.0, alternative_total);
                    let total: f64 = concentration.iter().sum();
                    let shared_constant = crate::genetics::lgamma(total + f64::from(copies))
                        - crate::genetics::lgamma(total);
                    let row =
                        mixed_row_for(copies, allele_count, 1.0, alternative_total, inbreeding);
                    let mass: f64 = row.iter().map(|p| (p.get() - shared_constant).exp()).sum();
                    assert!(
                        (mass - 1.0).abs() < 1e-9,
                        "ploidy {copies}, {allele_count} alleles, F {inbreeding}, alternative \
                         total {alternative_total}: the row carries {mass} rather than one unit \
                         of probability"
                    );
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
    /// gap is of order `1/Σα`. Measured here across four decades of total: each tenfold rise in the
    /// concentration shrinks the disagreement by about tenfold, so a tolerance that passed at
    /// `Σα = 1e6` for the wrong reason would not narrow with it.
    #[test]
    fn the_wright_agreement_closes_as_the_concentration_grows() {
        let frequency = 0.2;
        let inbreeding = 0.5;
        let (hom_reference, heterozygote, _) =
            crate::genetics::wright_genotype_log_priors(frequency, inbreeding);
        let wright_het_over_hom_ref = heterozygote - hom_reference;

        let mut previous_gap = f64::INFINITY;
        for total in [1e2, 1e3, 1e4, 1e5] {
            let alternative = total * frequency;
            let row = mixed_row_for(2, 2, total - alternative, alternative, inbreeding);
            let gap = ((row[1].get() - row[0].get()) - wright_het_over_hom_ref).abs();
            assert!(
                gap < previous_gap / 5.0,
                "at Σα {total} the gap to Wright was {gap}, not appreciably below the {previous_gap} \
                 at a tenth the concentration — the limit is not being approached"
            );
            previous_gap = gap;
        }
    }

    /// **The full-inbreeding limit** (spec §12 test 3). At `F = 1` every heterozygote is impossible
    /// and the two homozygotes stand in the ratio `α_ref : α_alt` — the identical-by-descent branch
    /// alone, which draws one allele and counts it twice.
    ///
    /// **The estimator never delivers this**, so it tests the mathematics at its edge rather than a
    /// case the caller meets: production clamps its inbreeding estimate at 0.99, and ng's newtype is
    /// to be tightened to `[0, 1)`. That is why it drives the bare-coefficient path — after the
    /// tightening, no `InbreedingF` will be able to carry a 1.
    ///
    /// The heterozygote lands at `−∞` rather than at a probability floor: `log(1 − 1)` is `−∞`, and
    /// nothing in this module floors a log-prior. After the loop rescales the row that is a weight
    /// of exactly zero, which is the wanted answer — an impossible genotype, not a very unlikely
    /// one.
    #[test]
    fn at_full_inbreeding_only_the_homozygotes_survive_and_they_stand_at_the_concentration_ratio() {
        for alternative_total in [1e-3, 1e-2, 0.25] {
            let row = mixed_row_for(2, 2, 1.0, alternative_total, 1.0);
            assert_eq!(
                row[1].get(),
                f64::NEG_INFINITY,
                "the heterozygote must be impossible at F = 1, got {}",
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
        fill_marginalized_log_priors(&mut row, inbreeding);

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
}
