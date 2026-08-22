//! The log-prior row: a handful of allele copies drawn from a population whose composition
//! is not known. That is the Dirichlet-multinomial — the locus's allele frequencies drawn
//! from a Dirichlet and averaged out, with this sample's copies drawn from what is left.
//!
//! **This file holds the random-mating half.** Every genotype gets it. The other half, the
//! branch for two copies that are the same one counted twice, is the inbreeding mixture that
//! plan step B2 wraps around it (`doc/devel/ng/spec/calling_priors.md` §3.2).

use crate::genetics::lgamma;
use crate::ng::calling::genotype_prior::PriorRow;
use crate::ng::types::LogProb;

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
}
